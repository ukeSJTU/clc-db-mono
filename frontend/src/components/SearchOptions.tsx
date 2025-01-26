import type React from "react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";

interface SearchOptionsProps {
  searchType: "smiles" | "file";
  setSearchType: (type: "smiles" | "file") => void;
  smiles: string;
  setSmiles: (smiles: string) => void;
  handleFileUpload: (event: React.ChangeEvent<HTMLInputElement>) => void;
  error: string;
}

export function SearchOptions({
  searchType,
  setSearchType,
  smiles,
  setSmiles,
  handleFileUpload,
  error,
}: SearchOptionsProps) {
  return (
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
    </div>
  );
}
