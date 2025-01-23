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
};

export default nextConfig;
