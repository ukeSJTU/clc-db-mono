from django.shortcuts import render
from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.views import APIView

from .utils import perform_clustering
from rdkit import Chem
from rdkit.Chem import AllChem
from e3fp.pipeline import fprints_from_sdf
import numpy as np
from faiss import IndexFlatIP
import os
import shutil
import time

class VectorSearchViewSet(viewsets.ViewSet):
    def get_morgan_fingerprint(self, smiles, radius=2, bits=1024):
        print(f"Generating Morgan fingerprint for SMILES: {smiles}")
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=bits)
        print(f"Generated fingerprint shape: {np.array(fp).shape}")
        return np.array(fp)

    def get_e3fp_fingerprint(self, sdf_file, bits=1024, radius_multiplier=2, rdkit_invariants=True):
        print(f"Generating E3FP fingerprint for file: {sdf_file}")
        fprint_params = {'bits': bits, 'radius_multiplier': radius_multiplier, 'rdkit_invariants': rdkit_invariants}
        fprint = fprints_from_sdf(sdf_file, fprint_params=fprint_params)
        print(f"Generated fingerprint shape: {np.array(fprint).shape}")
        vector = np.zeros(bits, dtype=np.int8)      
        vector[fprint[0].indices] = 1
        # return np.array(fprint)
        return vector


    @action(detail=False, methods=['post'])
    def search(self, request):
        print(f"Received search request with data: {request.data}")
        search_type = request.data.get('type')
        print(f"Search type: {search_type}")
        
        # Validate search type
        if search_type not in ['smiles', 'file']:
            return Response({'error': 'Invalid search type. Must be "smiles" or "file"'}, 
                          status=status.HTTP_400_BAD_REQUEST)
        
        if search_type == 'smiles':
            smiles = request.data.get('query')
            print(f"Processing SMILES query: {smiles}")
            if not smiles:
                return Response({'error': 'SMILES string is required'}, status=status.HTTP_400_BAD_REQUEST)
            
            try:
                query = self.get_morgan_fingerprint(smiles)
                print("Loading embeddings file...")
                embeddings = np.load('./data/embeddings_morgan.npy')
                print(f"Embeddings shape: {embeddings.shape}")
                
                index = IndexFlatIP(1024)
                print(1)
                index.add(embeddings)
                print("Added embeddings to index")
                
                distances, indices = index.search(query.reshape(1, -1), 10)
                print(f"Search results - indices: {indices}, distances: {distances}")
                
                results = [{'index': int(idx), 'distance': float(dist)} 
                          for idx, dist in zip(indices[0], distances[0])]
                
                return Response({'results': results}, status=status.HTTP_200_OK)
                
            except Exception as e:
                print(f"Error in SMILES search: {str(e)}")
                return Response({'error': str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
                
        elif search_type == 'file':
            file = request.FILES.get('file')
            print(f"Processing file upload: {file.name if file else 'No file'}")
            if not file:
                return Response({'error': 'File is required'}, status=status.HTTP_400_BAD_REQUEST)
            
            try:
                temp_dir = f'temp_{time.time()}'
                print(f"Creating temp directory: {temp_dir}")
                os.makedirs(temp_dir, exist_ok=True)
                file_path = os.path.join(temp_dir, file.name)
                
                with open(file_path, 'wb+') as destination:
                    for chunk in file.chunks():
                        destination.write(chunk)
                print(f"Saved uploaded file to: {file_path}")
                
                query = self.get_e3fp_fingerprint(file_path)
                print("Loading embeddings file...")
                embeddings = np.load('./data/embeddings_e3fp.npy')
                print(f"Embeddings shape: {embeddings.shape}")
                
                index = IndexFlatIP(1024)
                index.add(embeddings)
                print("Added embeddings to index")
                
                distances, indices = index.search(query.reshape(1, -1), 10)
                print(f"Search results - indices: {indices}, distances: {distances}")
                
                results = [{'index': int(idx), 'distance': float(dist)} 
                          for idx, dist in zip(indices[0], distances[0])]
                
                print(f"Cleaning up temp directory: {temp_dir}")
                shutil.rmtree(temp_dir)
                
                return Response({'results': results}, status=status.HTTP_200_OK)
                
            except Exception as e:
                print(f"Error in file search: {str(e)}")
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
                return Response({'error': str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
                
        else:
            print(f"Invalid search type: {search_type}")
            return Response({'error': 'Invalid search type'}, status=status.HTTP_400_BAD_REQUEST)

class SDFUploaderViewSet(viewsets.ViewSet):
    def create(self, request):
        uploaded_files = request.FILES.getlist("files")

        if not uploaded_files:
            return Response(
                {"error": "No files uploaded."}, status=status.HTTP_400_BAD_REQUEST
            )

        random_folder = str(time.time()).replace(".", "-")
        saved_folder = f"cluster/{random_folder}"
        os.makedirs(f"{saved_folder}", exist_ok=True)

        # Process and save the uploaded files
        saved_files = []
        for file in uploaded_files:
            # Save the file to a desired location
            file_path = f"cluster/{random_folder}/{file.name}"
            with open(file_path, "wb") as f:
                f.write(file.read())
            saved_files.append(file_path)

        # Log the received files
        print(f"Received files: {saved_files}")

        return Response(
            {
                "message": "Files uploaded successfully",
                "saved_folder": saved_folder,
                "files": saved_files,
            },
            status=status.HTTP_200_OK,
        )


class ClusteringViewSet(viewsets.ViewSet):
    @action(detail=False, methods=['post'])
    def cluster_by_category(self, request):
        category = request.data.get('category')
        if not category:
            return Response({'error': 'Category is required'}, status=status.HTTP_400_BAD_REQUEST)
            
        try:
            # Get SDF files for this category from database
            sdf_files = self.get_sdf_files_by_category(category)
            print(f"Found {len(sdf_files)} SDF files for category: {category}", sdf_files)
            if not sdf_files:
                return Response({'error': 'No SDF files found for this category'}, 
                              status=status.HTTP_404_NOT_FOUND)
                
            # Use default parameters from example clustering
            descriptor = "E3FP"
            bits = 1024
            radius = 1.5
            rdkit_inv = True
            reduction_method = "PCA"
            cluster_method = "K-Means"
            clusters = 5

            print("Clustering parameters:", descriptor, bits, radius, rdkit_inv, reduction_method, cluster_method, clusters)
            
            # Perform clustering
            result = perform_clustering(
                "data/random_category",
                descriptor,
                bits,
                radius,
                rdkit_inv,
                2,  # rdkit_radius
                False,  # rdkit_use_features
                False,  # rdkit_use_bond_types
                False,  # rdkit_use_chirality
                reduction_method,
                cluster_method,
                clusters,
                "lloyd",  # knn_algro
                0.25,  # eps
                5  # min_samples
            )

            print(f"Clustering complete for category: {category}")
            print(f"results: {result}")
            
            return Response({
                "message": f"Clustering performed for category: {category}",
                "results": result
            })
            
        except Exception as e:
            return Response({'error': str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    def get_sdf_files_by_category(self, category):
        """Helper method to get SDF files for a category from database"""
        # TODO: Implement actual database query
        category="random_category"
        # This is a placeholder implementation
        return [
            f'data/{category}/56-54-2.sdf',
            f'data/{category}/118-10-5.sdf',
            f'data/{category}/130-89-2.sdf',
            f'data/{category}/130-95-0.sdf',
            f'data/{category}/90-39-1.sdf'
        ]

    def create(self, request):
        saved_folder = request.data.get("saved_folder")
        descriptor = request.data.get("descriptor")
        bits = int(request.data.get("bits"))
        radius = float(request.data.get("radius"))
        rdkit_inv = request.data.get("rdkitInv") == "true"
        reduction_method = request.data.get("reductionMethod")
        cluster_method = request.data.get("clusterMethod")
        clusters = int(request.data.get("clusters"))
        knn_algro = request.data.get("knnAlgro")
        eps = float(request.data.get("eps"))
        min_samples = int(request.data.get("minSamples"))

        if not saved_folder:
            return Response(
                {"error": "Saved folder not provided."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        # Set default values for RDKit-related parameters if descriptor is "E3FP"
        if descriptor == "E3FP":
            rdkit_radius = 2
            rdkit_use_features = False
            rdkit_use_bond_types = False
            rdkit_use_chirality = False
        else:
            rdkit_radius = int(request.data.get("rdkitRadius", 2))
            rdkit_use_features = request.data.get("rdkitUseFeatures", "false") == "true"
            rdkit_use_bond_types = (
                request.data.get("rdkitUseBondTypes", "false") == "true"
            )
            rdkit_use_chirality = (
                request.data.get("rdkitUseChirality", "false") == "true"
            )

        # Perform clustering based on the selected options
        result = perform_clustering(
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
        )

        # Remove the saved folder to free up disk space
        try:
            shutil.rmtree(saved_folder)
        except OSError as e:
            print(f"Error: {saved_folder} : {e.strerror}")

        return Response(result, status=status.HTTP_200_OK)

    @action(detail=False, methods=["post"])
    def example(self, request):
        # Set the path to the predefined example SDF files
        example_folder = os.path.expanduser("~/Desktop/ChemNexus/cluster/")
        example_type = request.data.get("example_type", "default")

        # Set default parameters for the example
        descriptor = "E3FP"
        bits = 1024
        radius = 1.5
        rdkit_inv = True
        reduction_method = "PCA"
        cluster_method = "K-Means"
        clusters = 5
        knn_algro = "lloyd"
        eps = 0.25
        min_samples = 5

        if example_type == "default":
            descriptor = "E3FP"
            # TODO: check for missing params
            # ... (existing default parameters)
        elif example_type == "knn":
            descriptor = "RDKit"
            cluster_method = "K-Means"
            # ... (set other parameters for K-Means example)
        else:
            return Response(
                {"error": "Invalid example type."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        # Perform clustering using the example SDF files and default parameters
        result = perform_clustering(
            example_folder,
            descriptor,
            bits,
            radius,
            rdkit_inv,
            2,  # rdkit_radius
            False,  # rdkit_use_features
            False,  # rdkit_use_bond_types
            False,  # rdkit_use_chirality
            reduction_method,
            cluster_method,
            clusters,
            knn_algro,
            eps,
            min_samples,
        )

        return Response(result, status=status.HTTP_200_OK)
