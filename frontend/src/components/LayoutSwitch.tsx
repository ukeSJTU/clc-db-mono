import React, { memo } from "react";
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
    const newLayout = currentLayout === "grid" ? "table" : "grid";
    onToggleLayout(newLayout);
  };

  return (
    <div className="flex items-center gap-2" aria-label="Layout Switch">
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

export default memo(LayoutSwitch);
