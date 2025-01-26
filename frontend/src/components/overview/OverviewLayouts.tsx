import React, { memo } from "react";
import type { MoleculeProps } from "@/types/molecule";
import MoleculeTableRow from "./TableRow";
import MoleculeCard from "./Card";

/**
 * Props for both Table and Grid layouts
 */
interface LayoutProps {
  molecules: MoleculeProps[];
}

/**
 * Table layout
 */
export const MolecularTableLayout: React.FC<LayoutProps> = memo(
  ({ molecules }) => {
    return (
      <div className="shadow-md rounded-lg overflow-hidden border border-gray-200 dark:border-gray-800 bg-white dark:bg-gray-950">
        <table className="min-w-full divide-y divide-gray-200">
          <thead className="bg-gray-50 dark:bg-gray-700">
            <tr>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-600 dark:text-gray-300 uppercase">
                Name
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-600 dark:text-gray-300 uppercase">
                CAS ID
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-600 dark:text-gray-300 uppercase">
                Category
              </th>
              <th className="px-6 py-3 text-left text-xs font-medium text-gray-600 dark:text-gray-300 uppercase">
                Operations
              </th>
            </tr>
          </thead>
          <tbody className="bg-white dark:bg-gray-900 divide-y divide-gray-200 dark:divide-gray-700">
            {molecules.map((molecule) => (
              <MoleculeTableRow key={molecule.cas_id} molecule={molecule} />
            ))}
          </tbody>
        </table>
      </div>
    );
  }
);

MolecularTableLayout.displayName = "MolecularTableLayout";

/**
 * Grid layout
 */
export const MolecularGridLayout: React.FC<LayoutProps> = memo(
  ({ molecules }) => {
    return (
      <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 2xl:grid-cols-6 gap-6 mx-auto p-4">
        {molecules.map((molecule) => (
          <MoleculeCard key={molecule.cas_id} molecule={molecule} />
        ))}
      </div>
    );
  }
);

MolecularGridLayout.displayName = "MolecularGridLayout";
