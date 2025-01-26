import React, { memo } from "react";
import { useRouter } from "next/navigation";

import { MoleculeProps } from "@/types/molecule";

import { Button } from "@/components/ui/button";
import CategoryBadge from "@/components/CategoryBadge";
import { ZipDownloadButton } from "@/components/DownloadButtons";

/**
 * Table row sub-component for a single molecule
 */
interface MoleculeTableRowProps {
  molecule: MoleculeProps;
}

const MoleculeTableRow: React.FC<MoleculeTableRowProps> = ({ molecule }) => {
  const router = useRouter();

  const handleDetail = (casId: string) => {
    router.push(`/detail/${casId}`);
  };

  return (
    <tr className="hover:bg-gray-300/10">
      {/* Name column */}
      <td
        className="px-6 py-4 max-w-xs overflow-hidden text-ellipsis whitespace-nowrap"
        title={molecule.name}
      >
        {molecule.name}
      </td>
      {/* CAS ID column */}
      <td className="px-6 py-4 whitespace-nowrap">{molecule.cas_id}</td>
      {/* Categories column */}
      <td className="px-6 py-4 max-w-xs overflow-hidden text-ellipsis whitespace-nowrap">
        {molecule.category.length > 0 ? (
          molecule.category.map((cat) => (
            <CategoryBadge key={cat.id} category={cat} abbreviate={false} />
          ))
        ) : (
          <CategoryBadge category={{ name: "None" }} />
        )}
      </td>
      {/* Operations column */}
      <td className="px-6 py-4">
        <div className="flex flex-row justify-start space-x-2">
          <ZipDownloadButton
            molecules={[molecule]}
            sdfFiles={[
              `${process.env.NEXT_PUBLIC_STATIC}/all_sdfs/${molecule.cas_id}.sdf`,
            ]}
          />
          <Button
            variant="outline"
            onClick={() => handleDetail(molecule.cas_id)}
          >
            Detail
          </Button>
        </div>
      </td>
    </tr>
  );
};

export default memo(MoleculeTableRow);
