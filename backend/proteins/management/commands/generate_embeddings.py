import os
import numpy as np
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from e3fp.pipeline import fprints_from_sdf
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit import Chem
from rdkit.Chem import AllChem

class Command(BaseCommand):
    help = 'Generate embeddings for Morgan and e3fp from provided data and SDF files.'

    def add_arguments(self, parser):
        parser.add_argument(
            '-p', '--path', 
            type=str, 
            required=True, 
            help='Path to the CSV file (data.csv).'
        )
        parser.add_argument(
            '-s', '--sdf_dir', 
            type=str, 
            required=True, 
            help='Directory containing SDF files.'
        )
        parser.add_argument(
            '--reverse', 
            action='store_true', 
            help='Sort CAS IDs in descending order.'
        )

    def handle(self, *args, **options):
        csv_path = options['path']
        sdf_dir = options['sdf_dir']
        reverse_sort = options['reverse']

        # Check if CSV file exists
        if not os.path.isfile(csv_path):
            raise CommandError(f"CSV file not found: {csv_path}")

        # Check if SDF directory exists
        if not os.path.isdir(sdf_dir):
            raise CommandError(f"SDF directory not found: {sdf_dir}")

        # Load data
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            raise CommandError(f"Error loading CSV: {e}")

        if 'CAS' not in df.columns or 'SMILES' not in df.columns:
            raise CommandError("CSV must contain 'CAS' and 'SMILES' columns.")

        # Filter necessary columns
        df = df[['CAS', 'SMILES']]

        # Sort by CAS
        df = df.sort_values('CAS', ascending=not reverse_sort)
        
        # Generate Morgan embeddings
        morgan_embeddings = []
        for idx, row in df.iterrows():
            smiles = str(row['SMILES'])
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                fp_array = np.zeros((1,))
                ConvertToNumpyArray(fp, fp_array)
                morgan_embeddings.append(fp_array)

        morgan_embeddings = np.array(morgan_embeddings, dtype=np.float32)

        # Generate e3fp embeddings
        e3fp_embeddings = []
        for idx, row in df.iterrows():
            cas_id = row['CAS']
            sdf_file = os.path.join(sdf_dir, f"{cas_id}.sdf")
            if not os.path.isfile(sdf_file):
                self.stdout.write(f"Warning: SDF file not found for CAS {cas_id}. Skipping.")
                continue

            try:
                fprint_params = {'bits': 1024, 'radius_multiplier': 2, 'rdkit_invariants': True}
                fprint = fprints_from_sdf(sdf_file, fprint_params=fprint_params)
                vector = fprint[0].to_vector(sparse=False, dtype=float)
                e3fp_embeddings.append(vector)
            except Exception as e:
                self.stdout.write(f"Error processing SDF for CAS {cas_id}: {e}")
                continue

        e3fp_embeddings = np.array(e3fp_embeddings, dtype=np.float32)

        # Save embeddings
        os.makedirs('./backend/data', exist_ok=True)
        np.save('./backend/data/embeddings_morgan.npy', morgan_embeddings)
        np.save('./backend/data/embeddings_e3fp.npy', e3fp_embeddings)

        self.stdout.write(self.style.SUCCESS("Embeddings generated and saved successfully."))