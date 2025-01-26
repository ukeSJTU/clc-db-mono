"use client";

// This component is used to display a badge for a Category with dynamic colors and optional abbreviation.
// Clicking the badge navigates to the download page for the Category.

import React from "react";
import { Badge } from "@/components/ui/badge";
import { useRouter } from "next/navigation";

interface ClassTypeBadgeProps {
  category: { name: string }; // Expecting a single category object
  abbreviate?: boolean; // Optional: Control if the name should be abbreviated
  className?: string; // Optional: Accept additional CSS classes
}

const dynamicColor = (typeName: string): string => {
  const baseColors = [
    "bg-red-100 text-red-800",
    "bg-green-100 text-green-800",
    "bg-blue-100 text-blue-800",
    "bg-yellow-100 text-yellow-800",
  ];
  let hash = 0;
  for (let i = 0; i < typeName.length; i++) {
    hash = typeName.charCodeAt(i) + ((hash << 5) - hash);
  }
  return baseColors[Math.abs(hash) % baseColors.length];
};

const abbreviateName = (name: string): string => {
  return name
    .split(/\s+/)
    .map((word) => word[0])
    .join("")
    .toUpperCase();
};

const CategoryBadge: React.FC<ClassTypeBadgeProps> = ({
  category,
  abbreviate = true,
  className = "",
}) => {
  const router = useRouter();
  const colorClass = dynamicColor(category.name);
  const displayText = abbreviate
    ? abbreviateName(category.name)
    : category.name;

  const handleBadgeClick = () => {
    router.push(`/download/categories/${category.name}`);
  };

  return (
    <Badge
      className={`px-2 inline-flex text-xs font-semibold rounded-full cursor-pointer ${colorClass} text-nowrap ${className}`}
      title={category.name}
      onClick={handleBadgeClick}
      role="button" // Explicitly declare this as a button for accessibility
      tabIndex={0} // Make the badge focusable
      onKeyDown={(e) => {
        if (e.key === "Enter" || e.key === " ") {
          handleBadgeClick(); // Allow navigation via Enter or Space key
        }
      }}
    >
      {displayText}
    </Badge>
  );
};

export default CategoryBadge;
