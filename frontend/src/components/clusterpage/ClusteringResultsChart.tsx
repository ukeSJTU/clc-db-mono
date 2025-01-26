"use client";

import React, { useEffect, useRef } from "react";
import {
  Chart,
  ChartData,
  ChartOptions,
  registerables,
  TooltipItem,
  TooltipModel,
  TooltipPositionerFunction,
} from "chart.js";
import { Tooltip } from "chart.js";

declare module "chart.js" {
  interface TooltipPositionerMap {
    custom: TooltipPositionerFunction<"scatter">;
  }
}

export interface ClusteringResultsChartProps {
  results: {
    coordinates: string[][][];
    class_numbers: number[];
    ids: string[];
  };
}

// interface CustomTooltipItem extends TooltipItem<"scatter"> {
//   raw: {
//     x: number;
//     y: number;
//     cas_id: string;
//   };
// }

const ClusteringResultsChart: React.FC<ClusteringResultsChartProps> = ({
  results,
}) => {
  const chartRef = useRef<HTMLCanvasElement | null>(null);
  const chartInstanceRef = useRef<Chart | null>(null);

  useEffect(() => {
    if (chartRef.current) {
      const ctx = chartRef.current.getContext("2d");
      if (ctx) {
        Chart.register(...registerables);

        const data: ChartData<"scatter"> = {
          datasets: results.class_numbers.map((classNumber) => ({
            label: `Class ${classNumber}`,
            data: results.coordinates[classNumber].map((_, index) => ({
              x: parseFloat(results.coordinates[classNumber][index][0]),
              y: parseFloat(results.coordinates[classNumber][index][1]),
              cas_id: results.coordinates[classNumber][index][2], // Add cas_id here
            })),
            backgroundColor: `rgba(${Math.floor(
              Math.random() * 256
            )}, ${Math.floor(Math.random() * 256)}, ${Math.floor(
              Math.random() * 256
            )}, 0.7)`,
            pointRadius: 5,
            pointHoverRadius: 8,
          })),
        };

        // Register custom tooltip with proper typing
        Tooltip.positioners.custom = function (
          elements: readonly any[],
          position: { x: number; y: number }
        ) {
          return {
            x: position.x + 20,
            y: position.y - 20,
          };
        };

        const options: ChartOptions<"scatter"> = {
          responsive: true,
          plugins: {
            tooltip: {
              enabled: false,
              external: (context: {
                chart: Chart;
                tooltip: TooltipModel<"scatter"> & {
                  opacity: number;
                  dataPoints: { raw: { cas_id: string } }[];
                  caretX: number;
                  caretY: number;
                };
              }) => {
                const tooltipEl = document.getElementById("chartjs-tooltip");
                if (!tooltipEl) return;

                // Hide if no tooltip
                if (context.tooltip.opacity === 0) {
                  tooltipEl.style.opacity = "0";
                  return;
                }

                // Set tooltip content
                if (context.tooltip.dataPoints) {
                  const casId = context.tooltip.dataPoints[0].raw.cas_id;
                  tooltipEl.innerHTML = `
                            <div style="background: white; padding: 8px; border-radius: 4px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                              <div style="text-align: center; margin-bottom: 4px; font-weight: bold;">
                                ${casId}
                              </div>
                              <img src="https://compbio.sjtu.edu.cn/services/clc-db/static/2Dimages/${casId}.png" 
                                 alt="${casId}" 
                                 style="width: 200px; height: 200px; object-fit: contain; user-select: none;"/>
                            </div>
                          `;
                }

                // Position tooltip
                const { offsetLeft: positionX, offsetTop: positionY } =
                  context.chart.canvas;
                tooltipEl.style.opacity = "1";
                tooltipEl.style.left =
                  positionX + context.tooltip.caretX + "px";
                tooltipEl.style.top = positionY + context.tooltip.caretY + "px";
              },
              position: "custom",
            },
          },
          scales: {
            x: {
              type: "linear",
              grid: {
                display: true,
              },
            },
            y: {
              type: "linear",
              grid: {
                display: true,
              },
            },
          },
        };

        // Destroy the previous chart instance if it exists
        if (chartInstanceRef.current) {
          chartInstanceRef.current.destroy();
        }

        // Create a new chart instance
        chartInstanceRef.current = new Chart(ctx, {
          type: "scatter",
          data,
          options,
        });
      }
    }
  }, [results]);

  return (
    <div className="w-full h-96 relative">
      <canvas ref={chartRef} className="w-full h-full" />
      <div
        id="chartjs-tooltip"
        className="absolute opacity-0 pointer-events-none transition-opacity duration-200"
        style={{
          position: "absolute",
          pointerEvents: "none",
          transform: "translate(-50%, -100%)",
          userSelect: "none",
          WebkitUserSelect: "none",
        }}
      ></div>
    </div>
  );
};

export default ClusteringResultsChart;
