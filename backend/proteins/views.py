from django.shortcuts import render


from proteins.models import Category, Molecule
from django.db.models import Count
from rest_framework import filters, generics, pagination, viewsets, status
from rest_framework.decorators import action
from rest_framework.response import Response
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import pagination
from .filters import MoleculeFilter

from proteins.serializers import CategorySerializer, MoleculeSerializer


# custom pagination settings for overview api
class OverviewPagination(pagination.PageNumberPagination):
    page_size = 10
    page_size_query_param = "page_size"
    max_page_size = 100


class DownloadClassesPagination(pagination.PageNumberPagination):
    page_size = 10
    page_size_query_param = "page_size"
    max_page_size = 100


class DefaultPagination(pagination.PageNumberPagination):
    page_size = 10
    page_size_query_param = "page_size"
    max_page_size = 100


# Create your views here.
class OverviewViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    pagination_class = DefaultPagination


# api for searching with multiple query params single molecule data
class SearchViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer  # TODO: this serves search api, but should have its own serializer
    pagination_class = DefaultPagination

    filter_backends = (DjangoFilterBackend,)
    filterset_class = MoleculeFilter
