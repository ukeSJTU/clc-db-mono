from django.shortcuts import render

from rest_framework import viewsets, status
from rest_framework.decorators import action

# Create your views here.
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
