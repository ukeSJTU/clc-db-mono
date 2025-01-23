from django.core.management.base import BaseCommand
from proteins.models import Molecule, Category, Chirality
import pandas as pd
from django.db import transaction


class Command(BaseCommand):
    help = "Load molecule data into the database"
    fill_values = {
        "PubChem CID": 0,
        "SMILES": "N/A",
        "SMILES_IUPAC": "N/A",
        "Description": "N/A",
        "MF": "N/A",
        "MW": 0.0,
        "Heavy Atom Count": -1,
        "Ring Count": -1,
        "Hydrogen Bond Acceptor Count": -1,
        "Hydrogen Bond Donor Count": -1,
        "Rotatable Bond Count": -1,
        "types of chirality": "N/A",
        "InChI": "N/A",
        "InChIKey": "N/A",
        "Zero Point Correction": 0.0,
        "Thermal Correction Energy": 0.0,
        "Thermal Correction Enthalpy": 0.0,
        "Thermal Correction Gibbs": 0.0,
        "Sum Electronic Zero Point": 0.0,
        "Sum Electronic Thermal Energy": 0.0,
        "Sum Electronic Thermal Enthalpy": 0.0,
        "Sum Electronic Thermal Free Energy": 0.0,
        "HOMO Energy": 0.0,
        "LUMO Energy": 0.0,
        "HOMO-LUMO Gap": 0.0,
    }

    def add_arguments(self, parser):
        parser.add_argument(
            "-p",
            "--path",
            type=str,
            help="Path to the CSV file",
        )

    def handle(self, *args, **kwargs):
        csv_path = kwargs["path"]
        df = pd.read_csv(csv_path).fillna(self.fill_values)

        # Create or get categories and chirality types from the dataframe
        category_objs = {}
        for category_list in df["Category"].unique():
            categories = [cat.strip() for cat in category_list.split(",")]
            for category in categories:
                if category not in category_objs:
                    category_objs[category], _ = Category.objects.get_or_create(
                        name=category
                    )

        chirality_objs = {}
        for chirality in df["types of chirality"].unique():
            chirality_objs[chirality], _ = Chirality.objects.get_or_create(
                name=chirality
            )

        success_count = 0
        failure_count = 0
        failures = []

        for idx, row in df.iterrows():
            try:
                with transaction.atomic():
                    molecule, created = Molecule.objects.update_or_create(
                        cas_id=row["CAS"],
                        defaults={
                            "name": row["Name"],
                            "pubchem_cid": str(int(row.get("PubChem CID", 0))),
                            "url": row.get("URL", ""),
                            "pubchem_url": row.get("PubChemURL", ""),
                            "smiles": row.get("SMILES", ""),
                            "smiles_iupac": row.get("SMILES_IUPAC", ""),
                            "molecule_formula": row.get("MF", ""),
                            "molecular_weight": float(row.get("MW", 0)),
                            "description": row.get("Description", ""),
                            "InChI": row.get("InChI"),
                            "InChIKey": row.get("InChIKey"),
                            "heavy_atom_count": (
                                int(row.get("Heavy Atom Count", -1))
                                if row.get("Heavy Atom Count", -1) != -1
                                else None
                            ),
                            "ring_count": (
                                int(row.get("Ring Count", -1))
                                if row.get("Ring Count", -1) != -1
                                else None
                            ),
                            "hydrogen_bond_acceptor_count": (
                                int(row.get("Hydrogen Bond Acceptor Count", -1))
                                if row.get("Hydrogen Bond Acceptor Count", -1) != -1
                                else None
                            ),
                            "hydrogen_bond_donor_count": (
                                int(row.get("Hydrogen Bond Donor Count", -1))
                                if row.get("Hydrogen Bond Donor Count", -1) != -1
                                else None
                            ),
                            "rotatable_bond_count": (
                                int(row.get("Rotatable Bond Count", -1))
                                if row.get("Rotatable Bond Count", -1) != -1
                                else None
                            ),
                            "zero_point_correction": float(
                                row.get("Zero-point correction", 0)
                            ),
                            "thermal_correction_energy": float(
                                row.get("Thermal correction to Energy", 0)
                            ),
                            "thermal_correction_enthalpy": float(
                                row.get("Thermal correction to Enthalpy", 0)
                            ),
                            "thermal_correction_gibbs": float(
                                row.get("Thermal correction to Gibbs Free Energy", 0)
                            ),
                            "sum_electronic_zero_point": float(
                                row.get("Sum of electronic and zero-point Energies", 0)
                            ),
                            "sum_electronic_thermal_energy": float(
                                row.get("Sum of electronic and thermal Energies", 0)
                            ),
                            "sum_electronic_thermal_enthalpy": float(
                                row.get("Sum of electronic and thermal Enthalpies", 0)
                            ),
                            "sum_electronic_thermal_free_energy": float(
                                row.get(
                                    "Sum of electronic and thermal Free Energies", 0
                                )
                            ),
                            "homo_energy": float(row.get("HOMO Energy (eV)", 0)),
                            "lumo_energy": float(row.get("LUMO Energy (eV)", 0)),
                            "homo_lumo_gap": float(row.get("HOMO-LUMO Gap (eV)", 0)),
                        },
                    )

                    categories = [cat.strip() for cat in row["Category"].split(",")]
                    molecule.category.set(
                        [
                            category_objs[cat]
                            for cat in categories
                            if cat in category_objs
                        ]
                    )

                    chirality = row.get("types of chirality")
                    if chirality in chirality_objs:
                        molecule.chirality.set([chirality_objs[chirality]])

                    molecule.save()
                    success_count += 1
                    self.stdout.write(
                        self.style.SUCCESS(f"Successfully processed: {molecule.name}")
                    )

            except Exception as e:
                failure_count += 1
                failures.append((row["Name"], str(e)))
                self.stdout.write(
                    self.style.ERROR(f"Failed to process {row['Name']}: {str(e)}")
                )

        # Print summary
        self.stdout.write("\nImport Summary:")
        self.stdout.write(
            self.style.SUCCESS(f"Successfully processed: {success_count}")
        )
        self.stdout.write(self.style.ERROR(f"Failed to process: {failure_count}"))

        if failures:
            self.stdout.write("\nFailure Details:")
            for name, error in failures:
                self.stdout.write(f"- {name}: {error}")
