// This file handles the download functionality for the application.
// It provides functions to generate and download CSV and SDF files, as well as a combined ZIP file containing both.
// The `downloadFiles` function takes the type of file to download, an array of molecule objects, and an array of SDF file URLs as arguments.
// It then calls the appropriate function to generate and download the selected file type.
// Note that a "missing_files.txt" file is included in the ZIP archive if any SDF files are missing.

import { MoleculeProps } from "@/types/molecule";
import JSZip from "jszip";
import { saveAs } from "file-saver";

export type DownloadType = "csv" | "sdf" | "zip" | "mulliken" | "hirshfeld";

const escapeCsvField = (field: string): string => {
  if (field.includes(",") || field.includes('"') || field.includes("\n")) {
    return `"${field.replace(/"/g, '""')}"`;
  }
  return field;
};

const generateCsvContent = (molecules: MoleculeProps[]) => {
  const csvHeaders = [
    "Name",
    "CAS ID",
    "PubChem CID",
    "Category",
    "URL",
    "PubChem URL",
    "SMILES",
    "Chirality",
    "Description",
    "SMILES IUPAC",
    "Molecule Formula",
    "Molecular Weight",
    "Heavy Atom Count",
    "Ring Count",
    "Hydrogen Bond Acceptor Count",
    "Hydrogen Bond Donor Count",
    "Rotatable Bond Count",
    "Zero-point correction",
    "Thermal correction to Energy",
    "Thermal correction to Enthalpy",
    "Thermal correction to Gibbs Free Energy",
    "Sum of electronic and zero-point Energies",
    "Sum of electronic and thermal Energies",
    "Sum of electronic and thermal Enthalpies",
    "Sum of electronic and thermal Free Energies",
    "HOMO Energy (eV)",
    "LUMO Energy (eV)",
    "HOMO-LUMO Gap (eV)",
  ];

  const csvRows = molecules.map((molecule) => {
    const categories =
      molecule.category?.map((type) => type.name).join(", ") || "";
    const chiralities =
      molecule.chirality?.map((type) => type.name).join(", ") || "";
    return [
      escapeCsvField(molecule.name || ""),
      escapeCsvField(molecule.cas_id || ""),
      escapeCsvField(molecule.pubchem_cid || ""),
      escapeCsvField(categories),
      escapeCsvField(molecule.url || ""),
      escapeCsvField(molecule.pubchem_url || ""),
      escapeCsvField(molecule.smiles || ""),
      escapeCsvField(chiralities),
      escapeCsvField(molecule.description || ""),
      escapeCsvField(molecule.smiles_iupac || ""),
      escapeCsvField(molecule.molecule_formula || ""),
      molecule.molecular_weight?.toFixed(3) || "",
      molecule.heavy_atom_count?.toString() || "",
      molecule.ring_count?.toString() || "",
      molecule.hydrogen_bond_acceptor_count?.toString() || "",
      molecule.hydrogen_bond_donor_count?.toString() || "",
      molecule.rotatable_bond_count?.toString() || "",
      molecule.zero_point_correction?.toFixed(6) || "",
      molecule.thermal_correction_energy?.toFixed(6) || "",
      molecule.thermal_correction_enthalpy?.toFixed(6) || "",
      molecule.thermal_correction_gibbs?.toFixed(6) || "",
      molecule.sum_electronic_zero_point?.toFixed(6) || "",
      molecule.sum_electronic_thermal_energy?.toFixed(6) || "",
      molecule.sum_electronic_thermal_enthalpy?.toFixed(6) || "",
      molecule.sum_electronic_thermal_free_energy?.toFixed(6) || "",
      molecule.homo_energy?.toFixed(6) || "",
      molecule.lumo_energy?.toFixed(6) || "",
      molecule.homo_lumo_gap?.toFixed(6) || "",
    ];
  });

  return [csvHeaders, ...csvRows].map((row) => row.join(",")).join("\n");
};

const fetchSdfFiles = async (
  sdfFiles: string[],
  zip: JSZip
): Promise<string[]> => {
  const missingFiles: string[] = [];
  for (const [index, sdfUrl] of sdfFiles.entries()) {
    try {
      const response = await fetch(sdfUrl);
      if (!response.ok) {
        missingFiles.push(`Missing SDF: ${sdfUrl}`);
        continue;
      }
      const text = await response.text();
      // console.log(
      //   "Fetched length:",
      //   text.length,
      //   "content:",
      //   text.slice(0, 100)
      // );
      zip.file(`molecule_${index + 1}.sdf`, text);
    } catch {
      missingFiles.push(`Missing SDF: ${sdfUrl}`);
    }
  }
  return missingFiles;
};

const downloadFiles = async (
  type: DownloadType,
  molecules: MoleculeProps[],
  sdfFiles: string[]
): Promise<{ success: boolean; missingFiles: string[] }> => {
  const zip = new JSZip();
  let missingFiles: string[] = [];
  const defaultFileName =
    molecules.length === 1 ? molecules[0].cas_id : "molecules";

  try {
    if (type === "csv") {
      const csvContent = generateCsvContent(molecules);
      const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
      saveAs(blob, `${defaultFileName}.csv`);
    } else if (type === "sdf") {
      missingFiles = await fetchSdfFiles(sdfFiles, zip);
      const blob = await zip.generateAsync({ type: "blob" });
      saveAs(blob, `${defaultFileName}.zip`);
    } else if (type === "zip") {
      const csvContent = generateCsvContent(molecules);
      zip.file(`${defaultFileName}.csv`, csvContent);
      missingFiles = await fetchSdfFiles(sdfFiles, zip);
      const blob = await zip.generateAsync({ type: "blob" });
      saveAs(blob, `${defaultFileName}.zip`);
    } else if (type === "mulliken" || type === "hirshfeld") {
      for (const molecule of molecules) {
        const url = `${process.env.NEXT_PUBLIC_STATIC}/${type}/${molecule.cas_id}.chg`;
        try {
          const response = await fetch(url);
          if (!response.ok) {
            missingFiles.push(`Missing ${type} file: ${url}`);
            continue;
          }
          const text = await response.text();
          zip.file(`${molecule.cas_id}-${type}.chg`, text);
        } catch {
          missingFiles.push(`Missing ${type} file: ${url}`);
        }
      }
      const blob = await zip.generateAsync({ type: "blob" });
      saveAs(blob, `${defaultFileName}-${type}.zip`);
    } else {
      throw new Error("Unknown download type");
    }

    return { success: missingFiles.length === 0, missingFiles };
  } catch (error) {
    console.error("Download error:", error);
    return { success: false, missingFiles };
  }
};

export default downloadFiles;
