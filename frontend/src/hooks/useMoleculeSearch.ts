"use client";

import { useState, useEffect } from "react";
import { useDebouncedCallback } from "use-debounce";
import api from "@/utils/api";

import { MoleculeProps } from "@/types/molecule";

/** Defines a single search option in your UI (CAS, Name, etc.). */
export interface SearchOption {
  displayName: string;
  searchName: string; // "cas_id" | "name" | "smiles" ...
}

/** Pagination state structure. */
export interface PaginationState {
  page: number;
  pageSize: number;
  totalPages: number;
}

/** The shape of the data returned by the /search/molecules endpoint. */
interface SearchResponse {
  results: MoleculeProps[];
  count: number;
}

/**
 * A custom hook that encapsulates the logic for:
 * - Searching molecules by different fields (cas_id, name, smiles).
 * - Managing loading/error states.
 * - Handling pagination.
 * - Debouncing the search calls.
 */
export function useMoleculeSearch(
  initialSearchOption: string,
  initialPageSize = 12
) {
  const [searchOpt, setSearchOpt] = useState(initialSearchOption);
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<MoleculeProps[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState("");
  const [searchInitiated, setSearchInitiated] = useState(false);

  // Pagination state
  const [pagination, setPagination] = useState<PaginationState>({
    page: 1,
    pageSize: initialPageSize,
    totalPages: 0,
  });

  /**
   * Main search function, debounced to reduce server load
   * if the user types quickly.
   */
  const performSearch = useDebouncedCallback(
    async (
      searchQuery = query,
      searchOption = searchOpt,
      page = pagination.page,
      pageSize = pagination.pageSize
    ) => {
      if (!searchQuery.trim()) {
        // If query is empty, reset everything
        setResults([]);
        setPagination((prev) => ({ ...prev, totalPages: 0 }));
        return;
      }

      setIsLoading(true);
      setError("");

      try {
        const res = await api.get<SearchResponse>(
          `/search/molecules?${searchOption}=${searchQuery}&page=${page}&page_size=${pageSize}`
        );
        const { results: fetchedResults = [], count = 0 } = res.data;
        setResults(fetchedResults);

        setPagination((prev) => ({
          ...prev,
          totalPages: Math.ceil(count / pageSize),
        }));
      } catch (err: unknown) {
        setError("Failed to fetch molecules");
        console.error(err);
        setResults([]);
        setPagination((prev) => ({ ...prev, totalPages: 0 }));
      } finally {
        setIsLoading(false);
      }
    },
    300
  );

  /**
   * Call this when the user explicitly triggers a search
   * or when using examples, special searches, etc.
   */
  const handleSearch = (
    searchQuery = query,
    searchOption = searchOpt,
    page = 1,
    pageSize = pagination.pageSize
  ) => {
    setSearchInitiated(true);
    setPagination((prev) => ({ ...prev, page }));
    performSearch(searchQuery, searchOption, page, pageSize);
  };

  /**
   * Update pagination, then re-run the search for the new page/pageSize.
   */
  const updatePagination = (newPagination: Partial<PaginationState>) => {
    // Use functional update so we can see the old state as 'prev'
    setPagination((prev) => {
      // Merge 'newPagination' into the old state
      const updatedPagination = { ...prev, ...newPagination };

      // If a search has already been initiated, re-run the search
      if (searchInitiated) {
        performSearch(
          query,
          searchOpt,
          updatedPagination.page,
          updatedPagination.pageSize
        );
      }

      // Return the new pagination state
      return updatedPagination;
    });
  };

  // Re-run the search whenever the user changes page or pageSize
  // if a search was previously initiated
  useEffect(() => {
    if (searchInitiated) {
      performSearch();
    }
  }, [searchInitiated, performSearch]);

  return {
    searchOpt,
    setSearchOpt,
    query,
    setQuery,
    results,
    isLoading,
    error,
    pagination,
    updatePagination,
    handleSearch,
    setSearchInitiated,
  };
}
