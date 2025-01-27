"use client";

// This is the main navigation bar component that displays the main menu items.
import Link from "next/link";
import * as React from "react";

import { cn } from "@/lib/utils";
import { Icons } from "@/components/Icons";

import {
  NavigationMenu,
  NavigationMenuContent,
  NavigationMenuItem,
  NavigationMenuLink,
  NavigationMenuList,
  NavigationMenuTrigger,
  navigationMenuTriggerStyle,
} from "@/components/ui/navigation-menu";
import Image from "next/image";
import nextConfig from "@/utils/config";

const Navbar = () => {
  const basePath = nextConfig.basePath;
  return (
    <div className="flex justify-between flex-row w-full top-0 p-4 z-10 sticky bg-white">
      <NavigationMenu>
        <NavigationMenuList>
          <NavigationMenuItem>
            <NavigationMenuTrigger>Overview</NavigationMenuTrigger>
            <NavigationMenuContent>
              <ul className="grid gap-3 p-6 md:w-[400px] lg:w-[500px] lg:grid-cols-[.75fr_1fr]">
                <li className="row-span-3">
                  <NavigationMenuLink asChild>
                    <Link
                      href="/"
                      className="flex h-full w-full select-none flex-col justify-end rounded-md bg-gradient-to-b from-muted/50 to-muted p-6 no-underline outline-none focus:shadow-md"
                    >
                      <Icons.logo className="h-6 w-6" />
                      <div className="mb-2 mt-4 text-lg font-medium">
                        CLC-DB-Overview
                      </div>
                      <p className="text-sm leading-tight text-muted-foreground">
                        Overview of all the molecules in the database on CLC-DB.
                        Data rich, interactive and user-friendly.
                      </p>
                    </Link>
                  </NavigationMenuLink>
                </li>
                <ListItem href="/overview/card/1" title="Card">
                  View all the molecules in the database in a card format.
                </ListItem>
                <ListItem href="/overview/table/1" title="Table">
                  View the table of all the molecules in the database.
                </ListItem>
              </ul>
            </NavigationMenuContent>
          </NavigationMenuItem>

          <NavigationMenuItem>
            <Link href="/download/categories" passHref>
              <NavigationMenuLink className={navigationMenuTriggerStyle()}>
                Download
              </NavigationMenuLink>
            </Link>
          </NavigationMenuItem>

          <NavigationMenuItem>
            <Link href="/search" passHref>
              <NavigationMenuLink className={navigationMenuTriggerStyle()}>
                Search
              </NavigationMenuLink>
            </Link>
          </NavigationMenuItem>

          <NavigationMenuItem>
            <Link href="/cluster" passHref>
              <NavigationMenuLink className={navigationMenuTriggerStyle()}>
                Cluster
              </NavigationMenuLink>
            </Link>
          </NavigationMenuItem>

          <NavigationMenuItem>
            <Link href="/vector" passHref>
              <NavigationMenuLink className={navigationMenuTriggerStyle()}>
                Vector
              </NavigationMenuLink>
            </Link>
          </NavigationMenuItem>

          <NavigationMenuItem>
            <Link href="/statistics" passHref>
              <NavigationMenuLink className={navigationMenuTriggerStyle()}>
                Stats
              </NavigationMenuLink>
            </Link>
          </NavigationMenuItem>

          <NavigationMenuItem>
            <Link href="/help" passHref>
              <NavigationMenuLink className={navigationMenuTriggerStyle()}>
                Help
              </NavigationMenuLink>
            </Link>
          </NavigationMenuItem>

          <NavigationMenuItem>
            <Link href="/contact" passHref>
              <NavigationMenuLink className={navigationMenuTriggerStyle()}>
                Contact
              </NavigationMenuLink>
            </Link>
          </NavigationMenuItem>
        </NavigationMenuList>
      </NavigationMenu>

      <div className="text-lg font-bold flex items-center">
        <Link href="/" className="flex items-center">
          <Image
            src={`${basePath}/favicon.webp`}
            alt="Favicon"
            width={32}
            height={32}
            className="mr-2"
          />
          CLC-DB
        </Link>
      </div>
    </div>
  );
};

export default Navbar;

const ListItem = React.forwardRef<
  React.ElementRef<"a">,
  React.ComponentPropsWithoutRef<typeof Link> & { title: string }
>(({ className, title, children, ...props }, ref) => {
  return (
    <li>
      <NavigationMenuLink asChild>
        <Link
          ref={ref}
          className={cn(
            "block select-none space-y-1 rounded-md p-3 leading-none no-underline outline-none transition-colors hover:bg-accent hover:text-accent-foreground focus:bg-accent focus:text-accent-foreground",
            className
          )}
          {...props}
        >
          <div className="text-sm font-medium leading-none">{title}</div>
          <p className="line-clamp-2 text-sm leading-snug text-muted-foreground">
            {children}
          </p>
        </Link>
      </NavigationMenuLink>
    </li>
  );
});
ListItem.displayName = "ListItem";
