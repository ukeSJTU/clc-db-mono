/** @type {import('next').NextConfig} */
const nextConfig = {
  typescript: {
    ignoreBuildErrors: true,
  },
  images: {
    remotePatterns: [
      {
        protocol: "https",
        hostname: "pubchem.ncbi.nlm.nih.gov",
        port: "",
        pathname: "/image/**",
      },
      {
        protocol: "https",
        hostname: "compbio.sjtu.edu.cn",
      },
    ],
  },
  basePath: "/services/clc-db", // This is the base path of the website. For example, if you are deploying to https://www.example.com/website, then you would use "/website"
};

export default nextConfig;
