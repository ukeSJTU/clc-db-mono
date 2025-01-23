// This file provides a series of interactive charts that display the statistical info of the database
// The charts are created using the Chart.js library
"use client";

import { Bar, Pie } from "react-chartjs-2";
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend,
  ArcElement,
} from "chart.js";
import { useEffect, useState } from "react";
import api from "@/utils/api";

interface CategoryData {
  category__name: string;
  count: number;
}

interface ChiralityData {
  chirality__name: string;
  count: number;
}

interface WeightData {
  labels: string[];
  values: number[];
}

interface EnergyData {
  labels: string[];
  values: number[];
}

interface ChartData {
  labels: string[];
  datasets: {
    label?: string;
    data: number[];
    backgroundColor: string | string[];
  }[];
}

ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend,
  ArcElement
);

const CategoryDistChart = () => {
  const [chartData, setChartData] = useState<ChartData>({
    labels: [],
    datasets: [],
  });

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await api.get("/stats/category");
        const data = response.data as CategoryData[];

        console.log(data);

        const labels = data.map((item: CategoryData) => item.category__name);
        const counts = data.map((item: CategoryData) => item.count);

        setChartData({
          labels,
          datasets: [
            {
              label: "Number of Molecules",
              data: counts,
              backgroundColor: "rgba(75,192,192,0.6)",
            },
          ],
        });
      } catch (error) {
        console.error("Error fetching data:", error);
      }
    };

    fetchData();
  }, []);

  return (
    <div className="w-full h-full">
      <Bar
        data={chartData}
        options={{
          maintainAspectRatio: false,
          responsive: true,
          plugins: {
            legend: { position: "top" },
            title: {
              display: false,
              text: "Category Distribution",
            },
          },
        }}
      />
    </div>
  );
};

const ChiralityDistChart = () => {
  const [chartData, setChartData] = useState<ChartData>({
    labels: [],
    datasets: [],
  });

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await api.get("/stats/chirality");
        const data = response.data as ChiralityData[];

        console.log(data);

        const labels = data.map((item: ChiralityData) => item.chirality__name);
        const counts = data.map((item: ChiralityData) => item.count);

        setChartData({
          labels,
          datasets: [
            {
              data: counts,
              backgroundColor: [
                "#FF6384",
                "#36A2EB",
                "#FFCE56",
                "#8BC34A",
                "#E91E63",
              ],
            },
          ],
        });
      } catch (error) {
        console.error("Error fetching data:", error);
      }
    };

    fetchData();
  }, []);

  return (
    <div className="w-full h-full">
      <Pie
        data={chartData}
        options={{
          maintainAspectRatio: false,
          responsive: true,
          plugins: {
            legend: { position: "top" },
            title: {
              display: false,
              text: "Chirality Distribution",
            },
          },
        }}
      />
    </div>
  );
};

const WeightDistributionChart = () => {
  const [chartData, setChartData] = useState<ChartData>({
    labels: [],
    datasets: [],
  });

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await api.get("/stats/weights");
        const data = response.data as WeightData;
        setChartData({
          labels: data.labels,
          datasets: [
            {
              label: "Molecule Weight",
              data: data.values,
              backgroundColor: "rgba(75, 192, 192, 0.2)",
            },
          ],
        });
      } catch (error) {
        console.log("Error fetching weight distribution data", error);
      }
    };

    fetchData();
  }, []);

  return (
    <div className="w-full h-full">
      <Bar
        data={chartData}
        options={{
          maintainAspectRatio: false,
          responsive: true,
          plugins: {
            legend: { position: "top" },
            title: {
              display: false,
              text: "Molecule Weight Distribution",
            },
          },
        }}
      />
    </div>
  );
};

const HUMODistributionChart = () => {
  const [chartData, setChartData] = useState<ChartData>({
    labels: [],
    datasets: [],
  });

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await api.get("/stats/humo");
        const data = response.data as EnergyData;
        setChartData({
          labels: data.labels,
          datasets: [
            {
              label: "HUMO Energy",
              data: data.values,
              backgroundColor: "rgba(255, 99, 132, 0.2)",
            },
          ],
        });
      } catch (error) {
        console.log("Error fetching HUMO distribution data", error);
      }
    };

    fetchData();
  }, []);

  return (
    <div className="w-full h-full">
      <Bar
        data={chartData}
        options={{
          maintainAspectRatio: false,
          responsive: true,
          scales: {
            x: {
              ticks: {
                autoSkip: false,
                maxRotation: 90,
                minRotation: 90,
                callback: function (value: number | string, index: number) {
                  // Show only every 3rd label to reduce clutter
                  const numValue =
                    typeof value === "string" ? parseFloat(value) : value;
                  return index % 3 === 0 ? this.getLabelForValue(numValue) : "";
                },
              },
            },
          },
          plugins: {
            legend: { position: "top" },
            title: {
              display: false,
              text: "HUMO Energy Distribution",
            },
          },
        }}
      />
    </div>
  );
};

const LUMODistributionChart = () => {
  const [chartData, setChartData] = useState<ChartData>({
    labels: [],
    datasets: [],
  });

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await api.get("/stats/lumo");
        const data = response.data as EnergyData;
        setChartData({
          labels: data.labels,
          datasets: [
            {
              label: "LUMO Energy",
              data: data.values,
              backgroundColor: "rgba(54, 162, 235, 0.2)",
            },
          ],
        });
      } catch (error) {
        console.log("Error fetching LUMO distribution data", error);
      }
    };

    fetchData();
  }, []);

  return (
    <div className="w-full h-full">
      <Bar
        data={chartData}
        options={{
          maintainAspectRatio: false,
          responsive: true,
          scales: {
            x: {
              ticks: {
                autoSkip: false,
                maxRotation: 90,
                minRotation: 90,
                callback: function (value: number | string, index: number) {
                  // Show only every 3rd label to reduce clutter
                  const numValue =
                    typeof value === "string" ? parseFloat(value) : value;
                  return index % 3 === 0 ? this.getLabelForValue(numValue) : "";
                },
              },
            },
          },
          plugins: {
            legend: { position: "top" },
            title: {
              display: false,
              text: "LUMO Energy Distribution",
            },
          },
        }}
      />
    </div>
  );
};

export {
  ChiralityDistChart,
  CategoryDistChart,
  WeightDistributionChart,
  HUMODistributionChart,
  LUMODistributionChart,
};
