import React from "react";
import { Label } from "@/components/ui/label";
import { Switch } from "@/components/ui/switch";

interface LayoutSwitchProps {
  currentLayout: "grid" | "table";
  onToggleLayout: (layout: "grid" | "table") => void;
}

const LayoutSwitch: React.FC<LayoutSwitchProps> = ({
  currentLayout,
  onToggleLayout,
}) => {
  const handleToggle = () => {
    onToggleLayout(currentLayout === "grid" ? "table" : "grid");
  };

  return (
    <div
      className="flex justify-center items-center gap-2"
      aria-label="Layout Switch"
    >
      <Label htmlFor="layout-switch" className="cursor-pointer">
        Grid
      </Label>
      <Switch
        id="layout-switch"
        checked={currentLayout === "grid"}
        onCheckedChange={handleToggle}
        aria-checked={currentLayout === "grid"}
      />
      <Label htmlFor="layout-switch" className="cursor-pointer">
        Table
      </Label>
    </div>
  );
};

export default LayoutSwitch;
