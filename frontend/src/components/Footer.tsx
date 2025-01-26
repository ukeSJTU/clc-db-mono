import Link from "next/link";

const Footer = () => {
  return (
    <footer className="bg-gray-900 text-gray-400 py-12 px-4 md:px-6 lg:px-8">
      <div className="container mx-auto grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-8">
        {/* Overview Section */}
        <FooterSection title="Overview">
          <FooterLink href="/overview/card/1">Card</FooterLink>
          <FooterLink href="/overview/table/1">Table</FooterLink>
        </FooterSection>

        {/* Tools Section */}
        <FooterSection title="Tools">
          <FooterLink href="/download/categories">Download</FooterLink>
          <FooterLink href="/search">Search</FooterLink>
          <FooterLink href="/cluster">Cluster</FooterLink>
        </FooterSection>

        {/* Resources Section */}
        <FooterSection title="Resources">
          <FooterLink href="https://www.bidepharm.com/">毕得医药</FooterLink>
          <FooterLink href="https://pubchem.ncbi.nlm.nih.gov/">
            PubChem
          </FooterLink>
        </FooterSection>

        {/* Powered By Section */}
        <FooterSection title="Powered By">
          <FooterLink href="https://3dmol.org/doc/index.html" external>
            3Dmol.js
          </FooterLink>
          <FooterLink
            href="https://partridgejiang.github.io/Kekule.js/"
            external
          >
            Kekule.js
          </FooterLink>
          <FooterLink href="https://www.chartjs.org/" external>
            Chart.js
          </FooterLink>
        </FooterSection>
      </div>

      {/* Bottom Bar */}
      <div className="mt-8 text-center">
        <p>
          © 2024{" "}
          <Link
            href="https://compbio.sjtu.edu.cn/home.html"
            className="underline hover:text-gray-200 transition-colors"
          >
            Yang Lab
          </Link>
          . All rights reserved. The data is available under a CC-BY-NC 4.0
          license.
        </p>
      </div>
    </footer>
  );
};

export default Footer;

interface FooterSectionProps {
  title: string;
  children: React.ReactNode;
}

const FooterSection: React.FC<FooterSectionProps> = ({ title, children }) => (
  <div className="space-y-4">
    <h4 className="text-gray-300 font-semibold">{title}</h4>
    <ul className="space-y-2">{children}</ul>
  </div>
);

interface FooterLinkProps {
  href: string;
  children: React.ReactNode;
  external?: boolean;
}

const FooterLink: React.FC<FooterLinkProps> = ({
  href,
  children,
  external,
}) => (
  <li>
    <Link
      href={href}
      target={external ? "_blank" : undefined}
      rel={external ? "noopener noreferrer" : undefined}
      className="hover:text-gray-200 transition-colors"
    >
      {children}
    </Link>
  </li>
);
