import React from "react";
import { Slider } from "@/components/ui/slider";
import { Input } from "@/components/ui/input";

interface SliderWithValueProps {
  value: number;
  min: number;
  max: number;
  step: number;
  onChange: (value: number) => void;
}

const SliderWithValue: React.FC<SliderWithValueProps> = ({
  value,
  min,
  max,
  step,
  onChange,
}) => {
  const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = parseFloat(event.target.value);
    if (!isNaN(newValue) && newValue >= min && newValue <= max) {
      onChange(newValue);
    }
  };

  return (
    <div className="flex flex-col space-y-4">
      {/* Slider Row */}
      <div className="flex items-center space-x-4">
        {/* Min Value */}
        <div className="w-16 flex items-center justify-center bg-gray-100 rounded-md px-2 py-1">
          <span className="text-sm font-medium text-gray-600">{min}</span>
        </div>

        {/* Slider */}
        <div className="flex-1">
          <Slider
            value={[value]}
            min={min}
            max={max}
            step={step}
            onValueChange={(newValue) => onChange(newValue[0])}
            aria-label="slider"
          />
        </div>

        {/* Max Value */}
        <div className="w-16 flex items-center justify-center bg-gray-100 rounded-md px-2 py-1">
          <span className="text-sm font-medium text-gray-600">{max}</span>
        </div>
      </div>

      {/* Current Value Input */}
      <div className="flex justify-center">
        <div className="w-24">
          <Input
            type="number"
            value={value}
            onChange={handleInputChange}
            min={min}
            max={max}
            step={step}
            className="text-center"
          />
        </div>
      </div>
    </div>
  );
};

export default SliderWithValue;
