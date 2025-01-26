"use client";
import React, { useState } from "react";
import type { MoleculeProps } from "@/types/molecule";
import {
  MolecularGridLayout,
  MolecularTableLayout,
} from "@/components/overview/OverviewLayouts";
import {
  PaginationComponent,
  PaginationComponentProps,
} from "@/components/overview/Pagination";
import LayoutSwitch from "@/components/LayoutSwitch";

/**
 * Container Props:
 *
 * @param molecules – An array of molecule data
 * @param initialLayout – Either "grid" or "table", determines initial layout
 * @param useLayoutSwitch – Whether to show the switch to toggle layouts
 * @param paginationProps – Props for the pagination component
 * @param topLeftComponent – Any custom React node to display in the top-left area
 * @param children – Additional JSX to render at the bottom of the container
 */
interface OverviewContainerProps {
  molecules: MoleculeProps[];
  initialLayout?: "grid" | "table";
  useLayoutSwitch?: boolean;
  paginationProps: PaginationComponentProps;
  topLeftComponent?: React.ReactNode;
  children?: React.ReactNode;
}

const OverviewContainer: React.FC<OverviewContainerProps> = ({
  molecules,
  initialLayout = "grid",
  useLayoutSwitch = true,
  paginationProps,
  topLeftComponent,
  children,
}) => {
  const [layout, setLayout] = useState<"grid" | "table">(initialLayout);

  if (molecules.length === 0) {
    return (
      <div className="text-center text-gray-500 dark:text-gray-300 py-8">
        No results found.
      </div>
    );
  }

  return (
    <div className="bg-white dark:bg-gray-800 shadow-md rounded-lg p-4 sm:p-6">
      {/* Header: Top-left content and Layout Switch */}
      <div className="flex flex-col sm:flex-row justify-between items-center mb-6">
        <div className="mb-4 sm:mb-0">{topLeftComponent}</div>
        {useLayoutSwitch && (
          <LayoutSwitch currentLayout={layout} onToggleLayout={setLayout} />
        )}
      </div>

      {/* Main Content: Grid or Table Layout */}
      <div className="overflow-x-auto">
        {layout === "grid" ? (
          <MolecularGridLayout molecules={molecules} />
        ) : (
          <MolecularTableLayout molecules={molecules} />
        )}
      </div>

      {/* Pagination Controls */}
      <div className="mt-6 flex justify-center">
        <PaginationComponent {...paginationProps} />
      </div>

      {/* Additional optional content */}
      {children}
    </div>
  );
};

export default OverviewContainer;
