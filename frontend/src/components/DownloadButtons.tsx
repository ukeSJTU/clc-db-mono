"use client";

// This file contains multiple types of download buttons that can be used to download molecules in different formats.
// Note that, even if user choose to download .sdf file, they will get a .zip file containing the .sdf files.

import { useState } from "react";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { Button } from "@/components/ui/button";
import { useToast } from "@/hooks/use-toast";
import { ArchiveIcon, DownloadIcon, FileTextIcon } from "lucide-react";
import React from "react";
import downloadFiles, { DownloadType } from "@/lib/download";
import { MoleculeProps } from "@/types/molecule";

interface BaseDownloadButtonProps {
  type: DownloadType;
  molecules: MoleculeProps[];
  sdfFiles: string[];
  label?: string;
  icon?: React.ReactNode;
}

const BaseDownloadButton: React.FC<BaseDownloadButtonProps> = ({
  type,
  molecules,
  sdfFiles,
  label = "Download",
  icon = <DownloadIcon />,
}) => {
  const { toast } = useToast();

  const handleDownload = async () => {
    const { success, missingFiles } = await downloadFiles(
      type,
      molecules,
      sdfFiles
    );
    if (success) {
      toast({
        title: "Download Successful",
        description: `Your ${type.toUpperCase()} file(s) have been downloaded.`,
      });
    } else {
      toast({
        title: missingFiles.length > 0 ? "Missing Files" : "Download Failed",
        description:
          missingFiles.length > 0
            ? `The following files were missing: ${missingFiles.join(", ")}`
            : "An error occurred while downloading the file(s).",
        variant: missingFiles.length > 0 ? "warning" : "destructive",
      });
    }
  };

  return (
    <Button variant="outline" onClick={handleDownload} title={label}>
      <div className="flex items-center gap-2">
        {icon}
        <span>{label}</span>
      </div>
    </Button>
  );
};

// SelectDownloadButton is a download button that allows the user to select the download type
interface SelectDownloadButtonProps {
  molecules: MoleculeProps[];
  sdfFiles: string[];
}

const SelectDownloadButton: React.FC<SelectDownloadButtonProps> = ({
  molecules,
  sdfFiles,
}) => {
  const [downloadOption, setDownloadOption] = useState<DownloadType>("zip");
  const { toast } = useToast();

  const handleDownload = async () => {
    const { success, missingFiles } = await downloadFiles(
      downloadOption,
      molecules,
      sdfFiles
    );
    if (success) {
      toast({
        title: "Download Successful",
        description: `Your ${downloadOption.toUpperCase()} file(s) have been downloaded.`,
      });
    } else {
      toast({
        title: missingFiles.length > 0 ? "Missing Files" : "Download Failed",
        description:
          missingFiles.length > 0
            ? `The following files were missing: ${missingFiles.join(", ")}`
            : "An error occurred while downloading the file(s).",
        variant: missingFiles.length > 0 ? "warning" : "destructive",
      });
    }
  };

  return (
    <div className="flex items-center gap-2">
      <Select
        value={downloadOption}
        onValueChange={(value) => setDownloadOption(value as DownloadType)}
      >
        <SelectTrigger>
          <SelectValue placeholder="Download Both (ZIP)" />
        </SelectTrigger>
        <SelectContent>
          <SelectItem value="csv">CSV</SelectItem>
          <SelectItem value="sdf">SDF</SelectItem>
          <SelectItem value="mulliken">Mulliken</SelectItem>
          <SelectItem value="hirshfeld">Hirshfeld</SelectItem>
          <SelectItem value="zip">All</SelectItem>
        </SelectContent>
      </Select>
      <Button onClick={handleDownload} title="Download">
        <DownloadIcon />
      </Button>
    </div>
  );
};

// SingleMoleculeDownloadButton is a download button for a single molecule
interface SingleMoleculeDownloadButtonProps {
  molecule: MoleculeProps;
  sdfFile: string;
  type: "csv" | "sdf" | "zip";
}
const SingleMoleculeDownloadButton: React.FC<
  SingleMoleculeDownloadButtonProps
> = ({ molecule, sdfFile, type }) => {
  return (
    <BaseDownloadButton
      type={type}
      molecules={[molecule]}
      sdfFiles={[sdfFile]}
      label={`${type.toUpperCase()}`}
    />
  );
};

// BulkMoleculeDownloadButton is a download button for multiple molecules
interface BulkMoleculeDownloadButtonProps {
  molecules: MoleculeProps[];
  sdfFiles: string[];
  type: "csv" | "sdf" | "zip";
}
const BulkMoleculeDownloadButton: React.FC<BulkMoleculeDownloadButtonProps> = ({
  molecules,
  sdfFiles,
  type,
}) => {
  return (
    <BaseDownloadButton
      type={type}
      molecules={molecules}
      sdfFiles={sdfFiles}
      label={`${type.toUpperCase()}`}
    />
  );
};

// CsvDownloadButton is a download button for CSV files
interface CsvDownloadButtonProps {
  molecules: MoleculeProps[];
  sdfFiles: string[];
}
const CsvDownloadButton: React.FC<CsvDownloadButtonProps> = ({
  molecules,
  sdfFiles,
}) => {
  return (
    <BaseDownloadButton
      type="csv"
      molecules={molecules}
      sdfFiles={sdfFiles}
      label="CSV"
      icon={<FileTextIcon />}
    />
  );
};

// ZipDownloadButton is a download button for ZIP files
interface ZipDownloadButtonProps {
  molecules: MoleculeProps[];
  sdfFiles: string[];
}
const ZipDownloadButton: React.FC<ZipDownloadButtonProps> = ({
  molecules,
  sdfFiles,
}) => {
  return (
    <BaseDownloadButton
      type="zip"
      molecules={molecules}
      sdfFiles={sdfFiles}
      label="CSV+SDF"
      icon={<ArchiveIcon />}
    />
  );
};

export {
  SelectDownloadButton,
  BaseDownloadButton,
  SingleMoleculeDownloadButton,
  BulkMoleculeDownloadButton,
  CsvDownloadButton,
  ZipDownloadButton,
};
