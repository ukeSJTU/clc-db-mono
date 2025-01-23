from django.shortcuts import render
from django.db.models import Count

from rest_framework import viewsets
from rest_framework.response import Response
from collections import defaultdict
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


# class HUMODistributionViewSet(viewsets.ViewSet):
#     def get_humo_distribution(self):
#         # Define HUMO ranges with 1/3 intervals
#         humo_ranges = []
#         humo_bins = {}
        
#         # Generate ranges from -14 to -3 with 1/3 steps
#         start = -14.0
#         while start < -3:
#             end = round(start + 1/3, 2)
#             # Ensure range string format is correct
#             if start < end:
#                 range_str = f"{start:.2f}-{end:.2f}"
#                 humo_ranges.append(range_str)
#                 humo_bins[range_str] = 0
#                 start = end
#             else:
#                 break

#         # Classify molecules into HUMO bins
#         for molecule in Molecule.objects.all():
#             humo = molecule.homo_energy
#             if humo is not None:
#                 for range_str in humo_ranges:
#                     try:
#                         # Split and convert to float only if range_str is valid
#                         if '-' in range_str:
#                             start, end = map(float, range_str.split('-'))
#                             if start <= humo < end:
#                                 humo_bins[range_str] += 1
#                                 break
#                     except ValueError:
#                         continue

#         # Prepare data for the chart
#         data = {
#             "labels": humo_ranges,
#             "values": [humo_bins[range_] for range_ in humo_ranges],
#         }
#         return data

#     def list(self, request):
#         return Response(self.get_humo_distribution())


# class LUMODistributionViewSet(viewsets.ViewSet):
#     def get_lumo_distribution(self):
#         # Define LUMO ranges with 1/3 intervals
#         lumo_ranges = []
#         lumo_bins = {}
        
#         # Generate ranges from -5 to 3 with 1/3 steps
#         start = -5.0
#         while start < 3:
#             end = round(start + 1/3, 2)
#             # Ensure range string format is correct
#             if start < end:
#                 range_str = f"{start:.2f}-{end:.2f}"
#                 lumo_ranges.append(range_str)
#                 lumo_bins[range_str] = 0
#                 start = end
#             else:
#                 break

#         # Classify molecules into LUMO bins
#         for molecule in Molecule.objects.all():
#             lumo = molecule.lumo_energy
#             if lumo is not None:
#                 for range_str in lumo_ranges:
#                     try:
#                         # Split and convert to float only if range_str is valid
#                         if '-' in range_str:
#                             start, end = map(float, range_str.split('-'))
#                             if start <= lumo < end:
#                                 lumo_bins[range_str] += 1
#                                 break
#                     except ValueError:
#                         continue

#         # Prepare data for the chart
#         data = {
#             "labels": lumo_ranges,
#             "values": [lumo_bins[range_] for range_ in lumo_ranges],
#         }
#         return data

#     def list(self, request):
#         return Response(self.get_lumo_distribution())



class DistributionViewSet(viewsets.ViewSet):
    def generate_ranges(self, start, end, step):
        """Generate a list of range strings with specified step size."""
        ranges = []
        current = start
        while current < end:
            next_value = round(current + step, 2)
            if current < next_value:
                ranges.append(f"{current:.2f}~{next_value:.2f}")
            current = next_value
        return ranges

    def classify_data(self, ranges, data_list, data_field):
        """
        Classify data points into specified ranges.
        :param ranges: List of range strings (e.g., "start-end").
        :param data_list: Queryset or list of objects with numeric `data_field`.
        :param data_field: The field name to classify (e.g., 'lumo_energy').
        :return: Dictionary with counts per range.
        """
        bins = defaultdict(int)
        for obj in data_list:
            value = getattr(obj, data_field, None)
            if value is not None:
                for range_str in ranges:
                    low, high = map(float, range_str.split('~'))
                    if low <= value < high:
                        bins[range_str] += 1
                        break
        return bins

    def get_distribution(self, start, end, step, data_field):
        """
        Generate distribution data for a specific range and field.
        :param start: Start of the range.
        :param end: End of the range.
        :param step: Step size.
        :param data_field: The field to classify.
        :return: Dictionary with `labels` and `values`.
        """
        ranges = self.generate_ranges(start, end, step)
        bins = self.classify_data(ranges, Molecule.objects.all(), data_field)
        return {
            "labels": ranges,
            "values": [bins[range_] for range_ in ranges],
        }


class HUMODistributionViewSet(DistributionViewSet):
    def list(self, request):
        # Define range for HUMO distribution (-14 to -3)
        data = self.get_distribution(start=-14.0, end=-3.0, step=1/3, data_field="homo_energy")
        return Response(data)


class LUMODistributionViewSet(DistributionViewSet):
    def list(self, request):
        # Define range for LUMO distribution (-5 to 3)
        data = self.get_distribution(start=-5.0, end=3.0, step=1/3, data_field="lumo_energy")
        return Response(data)