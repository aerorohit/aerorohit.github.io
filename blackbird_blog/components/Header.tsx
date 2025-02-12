import Link from "next/link"

export default function Header() {
  return (
    <header className="border-b border-monokai-lighter py-4">
      <nav className="container mx-auto px-4 flex justify-end items-center">
        <Link href="/" className="mr-4 hover:text-monokai-yellow">
          Home
        </Link>
        <Link href="/about" className="hover:text-monokai-yellow">
          About
        </Link>
      </nav>
    </header>
  )
}

