from django.shortcuts import render
from django.db.models import Count

from rest_framework import viewsets
from rest_framework.response import Response
from proteins.models import Molecule, Category
from proteins.serializers import MoleculeSerializer, CategorySerializer

# Create your views here.


# api for statistics
class StatisticsViewSet(viewsets.ViewSet):
    def list(self, request):
        total_molecules = Molecule.objects.count()
        total_categories = Category.objects.count()

        return Response(
            {
                "total_molecules": total_molecules,
                "total_categories": total_categories,
            }
        )


# api for category object
class CategoryViewSet(viewsets.ViewSet):
    def list(self, request):
        categories = Category.objects.all()
        serializer = CategorySerializer(categories, many=True)
        return Response(serializer.data)

    def retrieve(self, request, pk=None):
        molecules = Molecule.objects.filter(class_type__id=pk)
        serialized_molecules = MoleculeSerializer(molecules, many=True)
        return Response(serialized_molecules.data)


# api for weight distribution stats
class WeightDistributionViewSet(viewsets.ViewSet):
    def weight_distribution(self):
        weight_ranges = ["0-100", "100-200", "200-300", "300-400", "400-500", "500+"]

        weight_bins = {
            "0-100": 0,
            "100-200": 0,
            "200-300": 0,
            "300-400": 0,
            "400-500": 0,
            "500+": 0,
        }

        # Adjust weights and classify into bins
        for molecule in Molecule.objects.all():
            weight = molecule.molecular_weight
            if weight is not None:
                if 0 < weight < 100:
                    weight_bins["0-100"] += 1
                elif weight < 200:
                    weight_bins["100-200"] += 1
                elif weight < 300:
                    weight_bins["200-300"] += 1
                elif weight < 400:
                    weight_bins["300-400"] += 1
                elif weight < 500:
                    weight_bins["400-500"] += 1
                else:
                    weight_bins["500+"] += 1

        # Prepare data for the chart
        data = {
            "labels": weight_ranges,
            "values": [weight_bins[range_] for range_ in weight_ranges],
        }

        return data

    def list(self, request):
        return Response(self.weight_distribution())


class ChiralityDistributionViewSet(viewsets.ViewSet):
    def get_chirality_distribution(self):
        distribution_data = list(
            Molecule.objects.values("chirality__name")
            .annotate(count=Count("id"))
            .order_by()
        )

        return distribution_data

    def list(self, request):
        return Response(self.get_chirality_distribution())


class CategoryDistributionViewSet(viewsets.ViewSet):
    def get_category_distribution(self):
        distribution_data = list(
            Molecule.objects.values("category__name")
            .annotate(count=Count("id"))
            .order_by()
        )
        return distribution_data

    def list(self, request):
        distribution_data = self.get_category_distribution()
        return Response(distribution_data)


class HUMODistributionViewSet(viewsets.ViewSet):
    def get_humo_distribution(self):
        # Define HUMO ranges
        humo_ranges = ["-10--5", "-5-0", "0-5", "5-10", "10-15", "15+"]

        humo_bins = {
            "-10--5": 0,
            "-5-0": 0,
            "0-5": 0,
            "5-10": 0,
            "10-15": 0,
            "15+": 0,
        }

        # Classify molecules into HUMO bins
        for molecule in Molecule.objects.all():
            humo = molecule.homo_energy
            if humo is not None:
                if -10 <= humo < -5:
                    humo_bins["-10--5"] += 1
                elif humo < 0:
                    humo_bins["-5-0"] += 1
                elif humo < 5:
                    humo_bins["0-5"] += 1
                elif humo < 10:
                    humo_bins["5-10"] += 1
                elif humo < 15:
                    humo_bins["10-15"] += 1
                else:
                    humo_bins["15+"] += 1

        # Prepare data for the chart
        data = {
            "labels": humo_ranges,
            "values": [humo_bins[range_] for range_ in humo_ranges],
        }
        return data

    def list(self, request):
        return Response(self.get_humo_distribution())


class LUMODistributionViewSet(viewsets.ViewSet):
    def get_lumo_distribution(self):
        # Define LUMO ranges
        lumo_ranges = ["-10--5", "-5-0", "0-5", "5-10", "10-15", "15+"]

        lumo_bins = {
            "-10--5": 0,
            "-5-0": 0,
            "0-5": 0,
            "5-10": 0,
            "10-15": 0,
            "15+": 0,
        }

        # Classify molecules into LUMO bins
        for molecule in Molecule.objects.all():
            lumo = molecule.lumo_energy
            if lumo is not None:
                if -10 <= lumo < -5:
                    lumo_bins["-10--5"] += 1
                elif lumo < 0:
                    lumo_bins["-5-0"] += 1
                elif lumo < 5:
                    lumo_bins["0-5"] += 1
                elif lumo < 10:
                    lumo_bins["5-10"] += 1
                elif lumo < 15:
                    lumo_bins["10-15"] += 1
                else:
                    lumo_bins["15+"] += 1

        # Prepare data for the chart
        data = {
            "labels": lumo_ranges,
            "values": [lumo_bins[range_] for range_ in lumo_ranges],
        }
        return data

    def list(self, request):
        return Response(self.get_lumo_distribution())
