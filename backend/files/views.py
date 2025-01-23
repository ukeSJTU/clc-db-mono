from django.shortcuts import render
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from proteins.models import Molecule, Category
from proteins.serializers import MoleculeSerializer
from proteins.views import DefaultPagination



# Create your views here.
class DownloadViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    pagination_class = DefaultPagination

    @action(detail=False, methods=["get"], url_path="all")
    def all_classes(self, request):
        queryset = Molecule.objects.all()

    @action(detail=False, methods=["get"], url_path="all")
    def all_classes(self, request):
        categories = Category.objects.all()
        paginated_categories = self.paginate_queryset(
            categories
        )  # Apply pagination to categories
        result = []

        if paginated_categories is not None:
            for category in paginated_categories:
                molecules = Molecule.objects.filter(class_type=category)
                serializer = self.get_serializer(molecules, many=True)
                result.append(
                    {"class_type": category.name, "molecules": serializer.data}
                )

            return self.get_paginated_response(result)  # Return paginated response
        else:
            # Handle the case if pagination is not applicable or error
            return Response({"error": "Pagination error or no categories found."})
