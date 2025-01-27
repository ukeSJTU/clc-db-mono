"use client";

import { useState, useEffect } from "react";
import api from "@/utils/api";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Button } from "@/components/ui/button";
import CategorySelector from "@/components/CategorySelector";
import ClusteringResultsChart from "@/components/clusterpage/ClusteringResultsChart";
import OverviewContainer from "@/components/OverviewContainer";
import { SearchOptions } from "@/components/SearchOptions";
import type { MoleculeProps } from "@/types/molecule";
import { Loader2 } from "lucide-react";

interface ClusterResult {
  coordinates: string[][][];
  class_numbers: number[];
  ids: string[];
}

interface SearchResult {
  molecule: MoleculeProps;
  distance: number;
}

export default function VectorPageClient() {
  const [activeTab, setActiveTab] = useState("faiss");
  const [searchType, setSearchType] = useState<"smiles" | "file">("smiles");
  const [smiles, setSmiles] = useState("");
  const [searchResult, setSearchResult] = useState<SearchResult[]>([]);
  const [currentPage, setCurrentPage] = useState(1);
  const itemsPerPage = 10;
  const [selectedCategory, setSelectedCategory] = useState("");
  const pageNumber = 0;
  const [paginationState, setPaginationState] = useState({
    page: pageNumber,
    pageSize: 12,
    totalPages: 0,
  });

  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState("");
  const [clusteringResults, setClusteringResults] =
    useState<ClusterResult | null>(null);
  const [isInitialLoading, setIsInitialLoading] = useState(true);

  useEffect(() => {
    const fetchInitialData = async () => {
      try {
        // Fetch any initial data you need here
        // For example, fetching categories for the CategorySelector
      } catch (error) {
        console.error("Error fetching initial data:", error);
      } finally {
        setIsInitialLoading(false);
      }
    };

    fetchInitialData();
  }, []);

  const handleSearch = async () => {
    if (!smiles && searchType === "smiles") {
      setError("Please enter a SMILES string");
      return;
    }

    setIsLoading(true);
    setError("");
    setSearchResult([]);

    try {
      const { data } = await api.post("/cluster/vector/search/search/", {
        type: searchType,
        query: smiles,
      });

      console.log(data);

      if (data.results && data.results.length > 0) {
        setSearchResult(data.results);
      } else {
        setSearchResult([]);
      }
    } catch (err) {
      if (err instanceof Error) {
        setError(err.message);
      } else {
        setError("An unknown error occurred");
      }
    } finally {
      setIsLoading(false);
    }
  };

  const handleFileUpload = async (
    event: React.ChangeEvent<HTMLInputElement>
  ) => {
    const file = event.target.files?.[0];
    if (!file) {
      setError("Please select a valid SDF file");
      return;
    }

    setIsLoading(true);
    setError("");
    setSearchResult([]);

    try {
      const formData = new FormData();
      formData.append("file", file);
      formData.append("type", "file");

      const { data } = await api.post(
        "/cluster/vector/search/search/",
        formData,
        {
          headers: {
            "Content-Type": "multipart/form-data",
          },
        }
      );

      if (data.results && data.results.length > 0) {
        setSearchResult(data.results);
      } else {
        setSearchResult([]);
      }
    } catch (err) {
      if (err instanceof Error) {
        setError(err.message);
      } else {
        setError("An unknown error occurred");
      }
    } finally {
      setIsLoading(false);
    }
  };

  const handleCluster = async () => {
    if (!selectedCategory) {
      setError("Please select a category");
      return;
    }

    setIsLoading(true);
    setError("");
    setClusteringResults(null);

    try {
      const { data } = await api.post(
        "/cluster/cluster_by_category/cluster_by_category/",
        {
          category: selectedCategory,
        }
      );

      if (data.results) {
        setClusteringResults(data.results);
      } else {
        setError("No clustering results found");
      }
    } catch (err) {
      if (err instanceof Error) {
        setError(err.message);
      } else {
        setError("An unknown error occurred");
      }
    } finally {
      setIsLoading(false);
    }
  };

  if (isInitialLoading) {
    return (
      <div className="flex items-center justify-center h-screen">
        <Loader2 className="h-8 w-8 animate-spin" />
      </div>
    );
  }

  return (
    <div className="mx-auto p-4">
      <h1 className="text-3xl font-bold mb-6 text-center">Vector Operations</h1>

      <Tabs
        value={activeTab}
        onValueChange={setActiveTab}
        className="max-w-6xl mx-auto"
      >
        <TabsList className="grid w-full grid-cols-2 gap-4 max-w-2xl mx-auto mb-6">
          <TabsTrigger value="faiss">FAISS Search</TabsTrigger>
          <TabsTrigger value="cluster">Cluster</TabsTrigger>
        </TabsList>

        <TabsContent value="faiss">
          <div className="space-y-6 max-w-3xl mx-auto">
            <SearchOptions
              searchType={searchType}
              setSearchType={setSearchType}
              smiles={smiles}
              setSmiles={setSmiles}
              handleFileUpload={handleFileUpload}
              error={error}
            />
            <div className="flex justify-center">
              <Button onClick={handleSearch} disabled={isLoading}>
                {isLoading ? "Searching..." : "Search"}
              </Button>
            </div>
            {error && (
              <div className="text-red-500 text-sm mt-2 text-center">
                {error}
              </div>
            )}
          </div>
        </TabsContent>

        <TabsContent value="cluster">
          <div className="space-y-6 max-w-3xl mx-auto">
            <CategorySelector
              selectedCategory={selectedCategory}
              onCategoryChange={setSelectedCategory}
            />
            <div className="flex justify-center">
              <Button onClick={handleCluster} disabled={isLoading}>
                {isLoading ? "Clustering..." : "Perform Clustering"}
              </Button>
            </div>
            {error && (
              <div className="text-red-500 text-sm mt-2 text-center">
                {error}
              </div>
            )}
          </div>

          {clusteringResults && (
            <div className="mt-8 bg-white p-6 rounded-lg shadow-md max-w-4xl mx-auto">
              <h2 className="text-2xl font-semibold mb-4 text-center">
                Clustering Results
              </h2>
              <ClusteringResultsChart results={clusteringResults} />
            </div>
          )}
        </TabsContent>
      </Tabs>
      {searchResult.length > 0 ? (
        <div className="mt-8 max-w-6xl mx-auto text-center">
          <h2 className="text-2xl font-semibold mb-4">Search Results</h2>
          <OverviewContainer
            molecules={searchResult.map((result) => result.molecule)}
            paginationProps={{
              page: currentPage,
              setPage: setCurrentPage,
              pageSize: itemsPerPage,
              setPageSize: (size) =>
                setPaginationState((prev) => ({ ...prev, pageSize: size })),
              totalPages: Math.ceil(searchResult.length / itemsPerPage),
            }}
          />
        </div>
      ) : (
        activeTab === "faiss" && (
          <p className="text-gray-500 mt-4 text-center">No results found</p>
        )
      )}
    </div>
  );
}
