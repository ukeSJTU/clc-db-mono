import type { Metadata } from "next";
import localFont from "next/font/local";
import "./globals.css";
import Navbar from "@/components/Navbar";
import Footer from "@/components/Footer";
import { Toaster } from "@/components/ui/toaster";

const geistSans = localFont({
  src: "./fonts/GeistVF.woff",
  variable: "--font-geist-sans",
  weight: "100 900",
});
const geistMono = localFont({
  src: "./fonts/GeistMonoVF.woff",
  variable: "--font-geist-mono",
  weight: "100 900",
});

export const metadata: Metadata = {
  title: "CLC-DB",
  description: "A collection of chemical resources for students and educators.",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <head>
        <link rel="icon" href="/services/clc-db/favicon.webp" />
        <link rel="preconnect" href="https://fonts.googleapis.com" />
        <link rel="preconnect" href="https://fonts.gstatic.com" />
        <link
          href="https://fonts.googleapis.com/css2?family=Andika:ital,wght@0,400;0,700;1,400;1,700&family=Inter:wght@100..900&display=swap"
          rel="stylesheet"
        />
      </head>
      <body
      // className={`antialiased`} // ${geistSans.variable} ${geistMono.variable}
      >
        <div className="flex flex-col min-h-screen justify-between">
          <div className="h-18">
            <Navbar />
          </div>
          <div className="top-18">
            <main className="mb-auto">{children}</main>
          </div>

          <Footer />
        </div>
        <Toaster />
      </body>
    </html>
  );
}
