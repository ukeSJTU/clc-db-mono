"use client";

import React, { useRef, useEffect, useState } from "react";
import { Skeleton } from "@/components/ui/skeleton";

type Molecule3DViewerProps = {
  casId: string; // CAS Registry Number as a prop
};

const Molecule3DViewer: React.FC<Molecule3DViewerProps> = ({ casId }) => {
  const viewerRef = useRef<HTMLDivElement | null>(null);
  const [data, setData] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(true);

  useEffect(() => {
    if (!viewerRef.current || !data) return;

    // Dynamically import 3Dmol.js to avoid loading it server-side
    const load3DViewer = async () => {
      try {
        const $3Dmol = await import("3dmol/build/3Dmol.js");
        const viewer = new $3Dmol.createViewer(viewerRef.current, {
          backgroundColor: "white",
        });

        viewer.addModel(data, "sdf");
        viewer.setStyle({}, { stick: {} });
        viewer.zoomTo();
        viewer.render();

        setIsLoading(false);
      } catch (e) {
        console.log(e);
        setError(
          "Failed to initialize the molecular viewer. See console for more details."
        );
        setIsLoading(false);
      }
    };

    load3DViewer();
  }, [data]);

  useEffect(() => {
    setIsLoading(true);

    const fetchData = async () => {
      try {
        const response = await fetch(
          `${process.env.NEXT_PUBLIC_STATIC}/all_sdfs/${casId}.sdf`,
          { mode: "no-cors" }
        );

        if (!response.ok) {
          throw new Error(
            response.status === 404
              ? `.sdf file for CAS ID ${casId} was not found.`
              : `Error fetching data: ${response.statusText}`
          );
        }

        const textData = await response.text();
        setData(textData);
        setError(null); // Clear any previous errors
      } catch (e) {
        setError(
          e instanceof Error ? e.message : "An unexpected error occurred."
        );
        setData(null);
      } finally {
        setIsLoading(false);
      }
    };

    fetchData();
  }, [casId]);

  return (
    <div
      style={{
        width: "100%",
        height: "100%",
        position: "relative",
      }}
      aria-label={`3D viewer for molecule with CAS ID ${casId}`}
    >
      {isLoading && <Skeleton className="h-full w-full" />}
      {!isLoading && error && (
        <div className="h-full w-full flex items-center justify-center bg-gray-100 text-gray-500">
          <p className="text-center">{error}</p>
        </div>
      )}
      {!isLoading && !error && (
        <div
          ref={viewerRef}
          style={{
            width: "100%",
            height: "100%",
          }}
        />
      )}
    </div>
  );
};

export default Molecule3DViewer;
