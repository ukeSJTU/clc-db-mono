import React, { memo } from "react";
import { useRouter } from "next/navigation";

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";

import { MoleculeProps } from "@/types/molecule";
import CategoryBadge from "@/components/CategoryBadge";
import MoleculeFormulaSpan from "@/components/MoleculeFormulaSpan";
import Molecule2DViewer from "@/components/Molecule2DViewer";
import { ZipDownloadButton } from "@/components/DownloadButtons";

interface MoleculeCardProps {
  molecule: MoleculeProps;
}

const MoleculeCard: React.FC<MoleculeCardProps> = ({ molecule }) => {
  const router = useRouter();
  const sdfFiles = [
    `${process.env.NEXT_PUBLIC_STATIC}/all_sdfs/${molecule.cas_id}.sdf`,
  ];

  const handleClick = () => {
    router.push(`/detail/${molecule.cas_id}`);
  };

  return (
    <div className="flex flex-col bg-white dark:bg-gray-800 shadow-md rounded-lg overflow-hidden w-full sm:w-64">
      <Card className="h-full flex flex-col">
        <CardHeader className="flex flex-col justify-between items-start p-4 border-b dark:border-gray-700">
          <div className="space-y-1">
            <CardTitle className="text-lg font-semibold line-clamp-4">
              {molecule.name}
            </CardTitle>
            <p className="text-sm text-gray-500">CAS ID: {molecule.cas_id}</p>
          </div>
          <div className="flex flex-wrap mt-2 space-y-1">
            {molecule.category?.map((cat) => (
              <CategoryBadge key={cat.id} category={cat} abbreviate={false} />
            ))}
          </div>
        </CardHeader>
        <CardContent className="p-4 flex-grow">
          <div className="grid grid-cols-2 gap-4">
            <div className="col-span-2 flex justify-between">
              <div className="space-y-2">
                <p className="text-sm text-gray-500">Formula</p>
                <MoleculeFormulaSpan formula={molecule.molecule_formula} />
              </div>
              <div className="space-y-2 text-right">
                <p className="text-sm text-gray-500">Weight</p>
                <p>
                  {typeof molecule.molecular_weight === "number"
                    ? molecule.molecular_weight.toFixed(3)
                    : "N/A"}
                </p>
              </div>
            </div>
            <div className="col-span-2 space-y-2">
              <p className="text-sm text-gray-500">2D Structure</p>
              <div className="w-full h-full relative flex justify-center">
                <Molecule2DViewer molecule={molecule} />
              </div>
            </div>
          </div>
        </CardContent>
        <div className="px-4 py-3 bg-gray-50 dark:bg-gray-700 text-right flex gap-2 justify-between">
          <ZipDownloadButton molecules={[molecule]} sdfFiles={sdfFiles} />
          <Button variant="outline" onClick={handleClick}>
            Detail
          </Button>
        </div>
      </Card>
    </div>
  );
};

export default memo(MoleculeCard);
