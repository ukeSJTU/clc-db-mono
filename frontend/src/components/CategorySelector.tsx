"use client";

import { Category } from "@/types/molecule";
import React from "react";
import useSWR from "swr";
import { Label } from "@/components/ui/label";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import CategoryBadge from "@/components/CategoryBadge";
import api from "@/utils/api";

interface CategorySelectorProps {
  selectedCategory: string;
  onCategoryChange: (category: string) => void;
}

// Fetcher function for useSWR
const fetchCategories = async (): Promise<Category[]> => {
  const { data } = await api.get("/categories/");
  return data;
};

const CategorySelector: React.FC<CategorySelectorProps> = ({
  selectedCategory,
  onCategoryChange,
}) => {
  // Fetch categories using SWR
  const {
    data: categories = [],
    error,
    isLoading,
  } = useSWR("/categories/", fetchCategories);

  if (isLoading) {
    return <div>Loading...</div>;
  }

  if (error) {
    return <div>Error loading categories. Please try again later.</div>;
  }

  return (
    <div className="space-x-4 space-y-2 flex flex-row items-center justify-center">
      <Label htmlFor="category" className="text-xl font-semibold text-nowrap">
        Select a category:
      </Label>
      <Select
        value={selectedCategory}
        onValueChange={(value) => onCategoryChange(value)}
      >
        <SelectTrigger>
          <SelectValue
            placeholder={
              <CategoryBadge
                category={{ name: "All categories" }}
                abbreviate={false}
              />
            }
          />
        </SelectTrigger>
        <SelectContent>
          <SelectItem value="all">
            <CategoryBadge
              category={{ name: "All categories" }}
              abbreviate={false}
            />
          </SelectItem>
          {categories.map((category) => (
            <SelectItem key={category.id} value={category.name}>
              <CategoryBadge category={category} abbreviate={false} />
            </SelectItem>
          ))}
        </SelectContent>
      </Select>
    </div>
  );
};

export default CategorySelector;
