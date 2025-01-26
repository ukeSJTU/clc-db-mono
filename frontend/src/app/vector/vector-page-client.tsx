"use client";

import { useState } from "react";
import api from "@/utils/api";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";
import CategorySelector from "@/components/CategorySelector";
import ClusteringResultsChart from "@/components/clusterpage/ClusteringResultsChart";

interface ClusterResult {
  labels: number[];
  coordinates: number[][];
  silhouetteScore?: number;
  clusterCenters?: number[][];
}

export default function VectorPageClient() {
  const [activeTab, setActiveTab] = useState("faiss");
  const [searchType, setSearchType] = useState<"smiles" | "file">("smiles");
  const [smiles, setSmiles] = useState("");
  const [searchResult, setSearchResult] = useState("");
  const [selectedCategory, setSelectedCategory] = useState("");

  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState("");
  const [clusteringResults, setClusteringResults] =
    useState<ClusterResult | null>(null);

  const handleSearch = async () => {
    if (!smiles) {
      setError("Please enter a SMILES string");
      return;
    }

    setIsLoading(true);
    setError("");
    setSearchResult("");

    try {
      const { data } = await api.post("/cluster/vector/search/search/", {
        type: "smiles",
        query: smiles,
      });

      if (data.results && data.results.length > 0) {
        const formattedResults = data.results
          .map(
            (result: any, index: number) =>
              `#${index + 1}: Index=${result.index}, Similarity=${(
                1 - result.distance
              ).toFixed(4)}`
          )
          .join("\n");

        setSearchResult(formattedResults);
      } else {
        setSearchResult("No results found");
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
    setSearchResult("");

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
        const formattedResults = data.results
          .map(
            (result: any, index: number) =>
              `#${index + 1}: Index=${result.index}, Similarity=${(
                1 - result.distance
              ).toFixed(4)}`
          )
          .join("\n");

        setSearchResult(formattedResults);
      } else {
        setSearchResult("No results found");
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

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-2xl font-bold mb-4">Vector Operations</h1>

      <Tabs value={activeTab} onValueChange={setActiveTab}>
        <TabsList className="grid w-full grid-cols-2">
          <TabsTrigger value="faiss">FAISS Search</TabsTrigger>
          <TabsTrigger value="cluster">Cluster</TabsTrigger>
        </TabsList>

        <TabsContent value="faiss">
          <div className="space-y-4">
            <div className="flex space-x-4 mb-4">
              <Button
                variant={searchType === "smiles" ? "default" : "outline"}
                onClick={() => setSearchType("smiles")}
              >
                Enter SMILES
              </Button>
              <Button
                variant={searchType === "file" ? "default" : "outline"}
                onClick={() => setSearchType("file")}
              >
                Upload SDF File
              </Button>
            </div>

            {searchType === "smiles" && (
              <div>
                <label
                  htmlFor="smiles"
                  className="block text-sm font-medium text-gray-700 mb-2"
                >
                  Enter SMILES String
                </label>
                <Input
                  id="smiles"
                  value={smiles}
                  onChange={(e) => setSmiles(e.target.value)}
                  placeholder="e.g. C1=CC=CC=C1"
                  className="w-full"
                />
                {error && searchType === "smiles" && !smiles && (
                  <p className="text-sm text-red-500 mt-2">
                    Please enter a valid SMILES string
                  </p>
                )}
              </div>
            )}

            {searchType === "file" && (
              <div>
                <label
                  htmlFor="sdf-upload"
                  className="block text-sm font-medium text-gray-700 mb-2"
                >
                  Upload SDF File
                </label>
                <div className="flex items-center space-x-4">
                  <Input
                    id="sdf-upload"
                    type="file"
                    onChange={handleFileUpload}
                    accept=".sdf"
                    className="w-full"
                  />
                </div>
                {error && searchType === "file" && (
                  <p className="text-sm text-red-500 mt-2">
                    Please select a valid SDF file
                  </p>
                )}
              </div>
            )}

            <Button onClick={handleSearch} disabled={isLoading}>
              {isLoading ? "Searching..." : "Search"}
            </Button>

            {error && <div className="text-red-500 text-sm mt-2">{error}</div>}

            <div>
              <label
                htmlFor="search-result"
                className="block text-sm font-medium text-gray-700"
              >
                Search Result
              </label>
              <Textarea
                id="search-result"
                value={searchResult}
                readOnly
                placeholder="Search results will appear here"
                className="h-32"
              />
            </div>
          </div>
        </TabsContent>

        <TabsContent value="cluster">
          <div className="space-y-4">
            <CategorySelector
              selectedCategory={selectedCategory}
              onCategoryChange={setSelectedCategory}
            />
            <Button onClick={handleCluster} disabled={isLoading}>
              {isLoading ? "Clustering..." : "Perform Clustering"}
            </Button>

            {error && <div className="text-red-500 text-sm mt-2">{error}</div>}

            {clusteringResults && (
              <div className="mt-8">
                <h2 className="text-xl font-semibold mb-4">
                  Clustering Results
                </h2>
                <ClusteringResultsChart results={clusteringResults} />
              </div>
            )}
          </div>
        </TabsContent>
      </Tabs>
    </div>
  );
}
