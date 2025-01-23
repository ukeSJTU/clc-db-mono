"""
URL configuration for clcDB project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""

from django.contrib import admin
from django.urls import path, include

from rest_framework.routers import DefaultRouter

from proteins.views import OverviewViewSet, SearchViewSet
from stats.views import (
    StatisticsViewSet,
    WeightDistributionViewSet,
    ChiralityDistributionViewSet,
    CategoryDistributionViewSet,
    CategoryViewSet,
    HUMODistributionViewSet,
    LUMODistributionViewSet,
)

router = DefaultRouter()

router.register(r"overview/card", OverviewViewSet, basename="overview-card")
router.register(r"overview/table", OverviewViewSet, basename="overview-table")

router.register(r"search/molecules", SearchViewSet, basename="search")

router.register(r"statistics", StatisticsViewSet, basename="statistics")
router.register(r"stats/weights", WeightDistributionViewSet, basename="stats-weights")
router.register(
    r"stats/chirality", ChiralityDistributionViewSet, basename="stats-chirality"
)
router.register(
    r"stats/category", CategoryDistributionViewSet, basename="stats-category"
)
router.register(r"categories", CategoryViewSet, basename="category")
router.register(r"stats/humo", HUMODistributionViewSet, basename="stats-humo")
router.register(r"stats/lumo", LUMODistributionViewSet, basename="stats-lumo")


urlpatterns = [path("admin/", admin.site.urls), path("api/", include(router.urls))]
