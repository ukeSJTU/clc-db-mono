"use client";
import { useState, useEffect } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import api from "@/utils/api";
import OverviewContainer from "@/components/OverviewContainer";
import CategorySelector from "@/components/CategorySelector";

const fetchMoleculeData = async (
  page: number,
  pageSize: number,
  category: string
) => {
  if (category === "all") {
    category = "";
  }
  const url = `/search/molecules/?page=${page}&page_size=${pageSize}${
    category ? `&category=${encodeURIComponent(category)}` : ""
  }`;
  console.log(url);
  const response = await api.get(url);
  const data = response.data;
  const totalPages = Math.ceil(data.count / pageSize);
  return { results: data.results, totalPages };
};

const IndexPage = ({ params }: { params: { pageNumber: string } }) => {
  const pageNumber = parseInt(params.pageNumber, 10);
  const searchParams = useSearchParams();
  const initialCategory = searchParams.get("category") || "";
  const [molecules, setMolecules] = useState([]);
  const [selectedCategory, setSelectedCategory] = useState(initialCategory);
  const [paginationState, setPaginationState] = useState({
    page: pageNumber,
    pageSize: 12,
    totalPages: 0,
  });

  const router = useRouter();

  useEffect(() => {
    const fetchData = async () => {
      const { results, totalPages } = await fetchMoleculeData(
        paginationState.page,
        paginationState.pageSize,
        selectedCategory
      );
      setMolecules(results);
      setPaginationState((prevState) => ({ ...prevState, totalPages }));
    };
    fetchData();
  }, [paginationState.page, paginationState.pageSize, selectedCategory]);

  const handleCategoryChange = (category: string) => {
    setSelectedCategory(category);
    setPaginationState((prevState) => ({ ...prevState, page: 1 }));
    router.push(`/overview/card/1?category=${encodeURIComponent(category)}`);
  };

  const handlePaginationChange = (
    newPaginationState: typeof paginationState
  ) => {
    setPaginationState(newPaginationState);
    const url = `/overview/card/${newPaginationState.page}${
      selectedCategory
        ? `?category=${encodeURIComponent(selectedCategory)}`
        : ""
    }`;
    router.push(url);
  };

  return (
    <div className="mx-auto p-4">
      <h2 className="text-3xl font-bold text-gray-800 dark:text-gray-200 text-center mb-8">
        Molecule Overview
      </h2>
      <OverviewContainer
        molecules={molecules}
        useLayoutSwitch={false}
        paginationProps={{
          page: paginationState.page,
          setPage: (page) =>
            handlePaginationChange({ ...paginationState, page }),
          pageSize: paginationState.pageSize,
          setPageSize: (pageSize) =>
            handlePaginationChange({
              ...paginationState,
              pageSize,
            }),
          totalPages: paginationState.totalPages,
        }}
        topLeftComponent={
          <CategorySelector
            selectedCategory={selectedCategory}
            onCategoryChange={handleCategoryChange}
          />
        }
      />
    </div>
  );
};

export default IndexPage;
