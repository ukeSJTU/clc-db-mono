"use client";

import { useState } from "react";
import api from "@/utils/api";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";

export default function VectorPageClient() {
  const [activeTab, setActiveTab] = useState("faiss");
  const [smiles, setSmiles] = useState("");
  const [searchResult, setSearchResult] = useState("");
  const [selectedCategory, setSelectedCategory] = useState("");

  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState("");

  const handleSearch = async () => {
    if (!smiles) {
      setError("Please enter a SMILES string");
      return;
    }

    setIsLoading(true);
    setError("");

    try {
      const { data } = await api.post("/cluster/vector/search/", {
        type: "smiles",
        query: smiles,
      });

      setSearchResult(JSON.stringify(data.results, null, 2));
    } catch (err) {
      if (err instanceof Error) {
        setError(err.message);
      } else {
        setError("An unknown error occurred");
      }
      setSearchResult("");
    } finally {
      setIsLoading(false);
    }
  };

  const handleFileUpload = async (
    event: React.ChangeEvent<HTMLInputElement>
  ) => {
    const file = event.target.files?.[0];
    if (!file) return;

    setIsLoading(true);
    setError("");

    try {
      const formData = new FormData();
      formData.append("file", file);

      const { data } = await api.post("/cluster/vector/search", formData, {
        headers: {
          "Content-Type": "multipart/form-data",
        },
      });

      setSearchResult(JSON.stringify(data.results, null, 2));
    } catch (err) {
      if (err instanceof Error) {
        setError(err.message);
      } else {
        setError("An unknown error occurred");
      }
      setSearchResult("");
    } finally {
      setIsLoading(false);
    }
  };

  const handleCluster = () => {
    console.log(`Performed clustering for category: ${selectedCategory}`);
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
            <div>
              <label
                htmlFor="smiles"
                className="block text-sm font-medium text-gray-700"
              >
                Enter SMILES
              </label>
              <Input
                id="smiles"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                placeholder="Enter SMILES string"
              />
            </div>

            <div>
              <label
                htmlFor="sdf-upload"
                className="block text-sm font-medium text-gray-700"
              >
                Upload SDF File
              </label>
              <Input
                id="sdf-upload"
                type="file"
                onChange={handleFileUpload}
                accept=".sdf"
              />
            </div>

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
            <div>
              <label
                htmlFor="category"
                className="block text-sm font-medium text-gray-700"
              >
                Select Category
              </label>
              <Select onValueChange={setSelectedCategory}>
                <SelectTrigger>
                  <SelectValue placeholder="Select a category" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="category1">Category 1</SelectItem>
                  <SelectItem value="category2">Category 2</SelectItem>
                  <SelectItem value="category3">Category 3</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <Button onClick={handleCluster}>Perform Clustering</Button>
          </div>
        </TabsContent>
      </Tabs>
    </div>
  );
}
