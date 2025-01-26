import React from "react";

interface MoleculeFormulaSpanProps {
  formula?: string; // Optional in case of empty or undefined input
}

const MoleculeFormulaSpan: React.FC<MoleculeFormulaSpanProps> = ({
  formula,
}) => {
  // Handle cases where the formula is undefined, empty, or invalid
  if (!formula?.trim() || formula === "N/A") {
    return <span>N/A</span>; // Show "N/A" for invalid input
  }

  // Regex to match chemical elements and their counts
  const formulaParts = Array.from(formula.matchAll(/([A-Z][a-z]*)(\d*)/g));

  return (
    <p>
      {formulaParts.map(([element, count], index) => (
        <span key={`${element}-${index}`}>
          {element}
          {count && <sub>{count}</sub>}
        </span>
      ))}
    </p>
  );
};

export default MoleculeFormulaSpan;
