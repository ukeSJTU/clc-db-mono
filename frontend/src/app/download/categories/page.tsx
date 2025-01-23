// This is the page to display all available class types for download

"use client";

import { Category, MoleculeProps } from "@/types/molecule";
import React, { useEffect, useState } from "react";
import api from "@/utils/api";
import CategoryBadge from "@/components/CategoryBadge";
import { Button } from "@/components/ui/button";
import downloadFiles from "@/lib/download";
import { useToast } from "@/hooks/use-toast";
import { Loader2 } from "lucide-react";

const CategoriesPage = () => {
  const [categories, setCategories] = useState<Category[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    // Fetch categories from your API
    const fetchCategories = async () => {
      const response = api.get("/categories/");
      const data = (await response).data;
      setCategories(data);
    };

    fetchCategories();
  }, []);

  const { toast } = useToast();

  const handleDownloadAll = async () => {
    setIsLoading(true);
    const allMolecules: MoleculeProps[] = [];
    const sdfFiles: string[] = [];
    let nextPage = 1;
    const allPageSize = 30;

    try {
      // Fetch all molecules across all categories
      let hasNext = true;
      while (hasNext) {
        const response = await api.get(
          `/search/molecules?page=${nextPage}&page_size=${allPageSize}`
        );
        const data = response.data;
        if (data.next === null) {
          hasNext = false;
        } else {
          allMolecules.push(...data.results);
          sdfFiles.push(
            ...data.results.map(
              (molecule: MoleculeProps) =>
                `${process.env.NEXT_PUBLIC_STATIC}/all_sdfs/${molecule.cas_id}.sdf`
            )
          );
          nextPage++;
        }
      }

      const { success, missingFiles } = await downloadFiles(
        "zip",
        allMolecules,
        sdfFiles
      );
      if (success) {
        toast({
          variant: "default",
          title: "Download Successful",
          description: "All molecules have been downloaded.",
        });
      } else if (missingFiles.length > 0) {
        toast({
          title: "Missing SDF Files",
          description: `The following SDF files were missing: ${missingFiles.join(
            ", "
          )}`,
        });
      } else {
        toast({
          variant: "destructive",
          title: "Download Failed",
          description: "An error occurred while downloading all molecules.",
        });
      }
    } catch (error) {
      toast({
        variant: "destructive",
        title: "Download Failed",
        description: "An error occurred while downloading all molecules.",
      });
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="flex flex-col items-center p-4">
      <h1 className="text-3xl font-bold mb-4 text-center">
        Available Class Types
      </h1>
      <p className="text-lg text-gray-600 mb-4 text-center">
        Click on a badge to view all the molecules under that category.
      </p>
      <div className="flex flex-wrap justify-center gap-4 max-w-4xl pb-4">
        {categories.map((category) => (
          <div
            key={category.id}
            className="transition duration-300 ease-in-out transform hover:-translate-y-1"
          >
            <CategoryBadge
              key={category.id}
              category={{ name: category.name }}
              abbreviate={false}
              className="text-base font-medium"
            />
          </div>
        ))}
      </div>
      <div className="w-full max-w-4xl mt-8">
        <Button
          onClick={handleDownloadAll}
          className="w-full py-6 text-lg font-bold bg-gradient-to-r from-blue-300 to-purple-300 hover:from-blue-400 hover:to-purple-400 text-gray-800 shadow-lg"
        >
          {isLoading ? <Loader2 className="mr-2 h-6 w-6 animate-spin" /> : null}
          Download All Molecules
        </Button>
      </div>
    </div>
  );
};

export default CategoriesPage;
