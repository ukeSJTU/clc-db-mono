"use client";

import React from "react";
import SearchBar from "@/components/searchpage/SearchBar";
import OverviewContainer from "@/components/OverviewContainer";
import SearchOptionsGroup from "@/components/searchpage/SearchOptionsGroup";
import { SearchHeading, SearchTip } from "@/components/searchpage/SearchText";
import SearchExamples from "@/components/searchpage/SearchExamples";

import {
  DrawStructureComponent,
  MultiCasIDSearchComponent,
} from "@/components/searchpage/SpecialSearch";

import { useMoleculeSearch } from "@/hooks/useMoleculeSearch";

/**
 * Available search fields for your UI: CAS ID, Name, SMILES
 * (plus any others you might add).
 */
const searchOptions = [
  { displayName: "CAS ID", searchName: "cas_id" },
  { displayName: "Name", searchName: "name" },
  { displayName: "SMILES", searchName: "smiles" },
];

export default function SearchPage() {
  // Our custom hook with default search option "cas_id" and default pageSize = 12
  const {
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
  } = useMoleculeSearch("cas_id", 12);

  // Example: If you want to automatically trigger a search once the user
  // has typed something (and they've setSearchInitiated to `true`),
  // you can rely on the debounced effect.
  // (We're leaving it manual for now, so we only search on button or examples.)

  /**
   * A helper for special or example searches
   */
  const handleSpecialSearchInput = (
    queryStr: string,
    searchOptValue: string,
    triggerSearch = true
  ) => {
    setQuery(queryStr);
    setSearchOpt(searchOptValue);
    if (triggerSearch) {
      handleSearch(queryStr, searchOptValue, 1, pagination.pageSize);
    }
  };

  const handleExampleClick = (exampleQuery: string, exampleOpt: string) => {
    handleSpecialSearchInput(exampleQuery, exampleOpt, true);
  };

  // If you prefer a separate error UI, you can do so here
  // or show it inline near the search bar

  return (
    <div className="flex flex-col items-center py-12 space-y-4">
      {/* Heading & Intro */}
      <div className="flex flex-col gap-2 w-full max-w-md sm:max-w-lg md:max-w-2xl">
        <SearchHeading />
        <SearchTip />

        {/* Search Option Checkboxes (CAS, Name, SMILES) */}
        <SearchOptionsGroup
          options={searchOptions}
          setSearchOpt={setSearchOpt}
          searchOpt={searchOpt}
        />

        {/* Main search bar */}
        <SearchBar
          query={query}
          setQuery={setQuery}
          handleSearch={() =>
            handleSearch(query, searchOpt, 1, pagination.pageSize)
          }
          onSmilesInput={(smiles) => {
            setQuery(smiles);
            setSearchOpt("smiles");
          }}
        />

        {/* Common Examples */}
        <SearchExamples onExampleClick={handleExampleClick} />

        {/* Additional specialized search UI:
            Kekule drawing + Multi-CAS file upload
         */}
        <div className="flex flex-row items-center justify-center px-4">
          <div className="flex justify-end mr-2 w-1/2">
            <DrawStructureComponent
              onSubmit={(smiles) => handleSpecialSearchInput(smiles, "smiles")}
              onClose={() => setSearchInitiated(true)}
            />
          </div>
          <div className="flex justify-start ml-2 w-1/2">
            <MultiCasIDSearchComponent
              onSubmit={(casIds) => handleSpecialSearchInput(casIds, "cas_id")}
              onClose={() => setSearchInitiated(true)}
            />
          </div>
        </div>
      </div>

      {/* If there's an error from the search hook, show it */}
      {error && <div className="text-red-500 text-sm mt-2">{error}</div>}

      {/* Display loading state if needed */}
      {isLoading && (
        <div className="flex justify-center mt-2">
          <p>Loading...</p>
        </div>
      )}

      {/* Results & Pagination */}
      <OverviewContainer
        molecules={results}
        paginationProps={{
          page: pagination.page,
          setPage: (page) => updatePagination({ page }),
          pageSize: pagination.pageSize,
          setPageSize: (pageSize) => updatePagination({ pageSize, page: 1 }),
          totalPages: pagination.totalPages,
        }}
      />
    </div>
  );
}
