import os
import numpy as np
from e3fp.pipeline import fprints_from_sdf, mol_from_sdf
from rdkit import Chem
from rdkit.Chem import AllChem
from .embedding import embedding, cluster


def perform_clustering(
    saved_folder,
    descriptor,
    bits,
    radius,
    rdkit_inv,
    rdkit_radius,
    rdkit_use_features,
    rdkit_use_bond_types,
    rdkit_use_chirality,
    reduction_method,
    cluster_method,
    clusters,
    knn_algro,
    eps,
    min_samples,
):
    print(f"Starting clustering with folder: {saved_folder}")
    fprint_list = []
    id_list = []
    
    print(f"Processing files in directory: {os.listdir(saved_folder)}")
    
    for file_name in os.listdir(saved_folder):
        file_path = os.path.join(saved_folder, file_name)
        print(f"Processing file: {file_path}")
        
        if os.path.isfile(file_path) and file_name.endswith(".sdf"):
            try:
                mol = mol_from_sdf(file_path)
                print(f"Successfully loaded molecule from: {file_name}")
                
                if descriptor == "E3FP":
                    print("Generating E3FP fingerprint")
                    fprint = fprints_from_sdf(
                        file_path,
                        fprint_params={
                            "bits": bits,
                            "radius_multiplier": radius,
                            "rdkit_invariants": rdkit_inv,
                        },
                    )
                    vector = fprint[0].to_vector(sparse=False, dtype=int)
                elif descriptor == "RDKit":
                    print("Generating RDKit fingerprint")
                    vector = AllChem.GetMorganFingerprintAsBitVect(
                        mol,
                        radius=rdkit_radius,
                        nBits=bits,
                        useFeatures=rdkit_use_features,
                        useBondTypes=rdkit_use_bond_types,
                        useChirality=rdkit_use_chirality,
                    )
                fprint_list.append(vector)
                id = os.path.splitext(file_name)[0]
                id_list.append(id)
            except AttributeError as e:
                print(f"AttributeError processing {file_name}: {str(e)}")
                continue
            except Exception as e:
                print(f"Unexpected error processing {file_name}: {str(e)}")
                continue

    print(f"Total fingerprints generated: {len(fprint_list)}")
    
    try:
        print("Starting embedding process")
        features = embedding(np.array(fprint_list), reduction_method)
        print(f"Embedding complete. Shape: {features.shape}")
        
        print("Starting clustering process")
        coords_with_id_list, class_num = cluster(
            features,
            id_list,
            cluster_method,
            {
                "n_cluster": clusters,
                "knn_algro": knn_algro,
                "eps": eps,
                "min_samples": min_samples,
            },
        )
        print(f"Clustering complete. Number of clusters: {len(set(class_num))}")
    except Exception as e:
        print(f"Error during embedding or clustering: {str(e)}")
        raise

    # Prepare the clustering results as JSON data
    result = {
        "coordinates": [coords.tolist() for coords in coords_with_id_list],
        "class_numbers": list(class_num),
        "ids": id_list,
    }
    print("Successfully created result dictionary")
    print(f"Result dictionary: {result}")
    return result
