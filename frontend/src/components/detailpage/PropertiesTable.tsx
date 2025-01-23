import React from "react";
import { MoleculeProps } from "@/types/molecule";
import MoleculeFormulaSpan from "@/components/MoleculeFormulaSpan";

const MoleculePropertiesTable = ({ molecule }: { molecule: MoleculeProps }) => {
  const formatCount = (count: number | null | undefined) => {
    if (count === null || count === undefined || count === -1) {
      return "N/A";
    } else {
      return count;
    }
  };

  return (
    <div className="overflow-x-auto">
      <div className="mb-4">
        <h3 className="text-xl font-semibold">Properties</h3>
        <small className="text-xs text-gray-500 dark:text-gray-400">
          Note: the properties below are calculated by rdkit.
        </small>
      </div>
      <table className="w-full table-auto">
        <tbody>
          <tr>
            <td className="py-2 pr-4 font-semibold">IUPAC Name</td>
            <td className="py-2">{molecule.smiles_iupac || "N/A"}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">InChI</td>
            <td className="py-2">{molecule.InChI || "N/A"}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">InChIKey</td>
            <td className="py-2">{molecule.InChIKey || "N/A"}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 fon</small>t-semibold">
              Molecular Formula
            </td>
            <td className="py-2">
              <MoleculeFormulaSpan formula={molecule.molecule_formula} />
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">Molecular Weight</td>
            <td className="py-2">{molecule.molecular_weight || "N/A"}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">Heavy Atom Count</td>
            <td className="py-2">{formatCount(molecule.heavy_atom_count)}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">Ring Count</td>
            <td className="py-2">{formatCount(molecule.ring_count)}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">H-Bond Acceptor Count</td>
            <td className="py-2">
              {formatCount(molecule.hydrogen_bond_acceptor_count)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">H-Bond Donor Count</td>
            <td className="py-2">
              {formatCount(molecule.hydrogen_bond_donor_count)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">Rotatable Bond Count</td>
            <td className="py-2">
              {formatCount(molecule.rotatable_bond_count)}
            </td>
          </tr>

          <div className="mb-4">
            <small className="text-xs text-gray-500 dark:text-gray-400 block mt-1">
              Note: the following properties are calculated by DFT.
            </small>
          </div>
          <tr>
            <td className="py-2 pr-4 font-semibold">Zero Point Correction</td>
            <td className="py-2">
              {formatCount(molecule.zero_point_correction)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">
              Thermal Correction Energy
            </td>
            <td className="py-2">
              {formatCount(molecule.thermal_correction_energy)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">
              Thermal Correction Enthalpy
            </td>
            <td className="py-2">
              {formatCount(molecule.thermal_correction_enthalpy)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">
              Thermal Correction Gibbs
            </td>
            <td className="py-2">
              {formatCount(molecule.thermal_correction_gibbs)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">
              Sum Electronic + Zero Point
            </td>
            <td className="py-2">
              {formatCount(molecule.sum_electronic_zero_point)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">
              Sum Electronic + Thermal Energy
            </td>
            <td className="py-2">
              {formatCount(molecule.sum_electronic_thermal_energy)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">
              Sum Electronic + Thermal Enthalpy
            </td>
            <td className="py-2">
              {formatCount(molecule.sum_electronic_thermal_enthalpy)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">
              Sum Electronic + Thermal Free Energy
            </td>
            <td className="py-2">
              {formatCount(molecule.sum_electronic_thermal_free_energy)}
            </td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">HOMO Energy</td>
            <td className="py-2">{formatCount(molecule.homo_energy)}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">LUMO Energy</td>
            <td className="py-2">{formatCount(molecule.lumo_energy)}</td>
          </tr>
          <tr>
            <td className="py-2 pr-4 font-semibold">HOMO-LUMO Gap</td>
            <td className="py-2">{formatCount(molecule.homo_lumo_gap)}</td>
          </tr>
        </tbody>
      </table>
    </div>
  );
};

export default MoleculePropertiesTable;
