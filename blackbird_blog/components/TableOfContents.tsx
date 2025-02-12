"use client"

import { useEffect, useState } from "react"
import { cn } from "@/lib/utils"

interface TableOfContentsProps {
  items: Array<{
    id: string
    title: string
    level: number
  }>
}

export default function TableOfContents({ items }: TableOfContentsProps) {
  const [activeId, setActiveId] = useState<string>("")

  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            setActiveId(entry.target.id)
          }
        })
      },
      { rootMargin: "0% 0% -80% 0%" },
    )

    items.forEach(({ id }) => {
      const element = document.getElementById(id)
      if (element) {
        observer.observe(element)
      }
    })

    return () => observer.disconnect()
  }, [items])

  return (
    <nav className="sticky top-4">
      <h2 className="text-lg font-semibold mb-4 text-monokai-orange">Contents</h2>
      <ul className="space-y-2">
        {items.map(({ id, title, level }) => (
          <li key={id} style={{ paddingLeft: `${(level - 1) * 1}rem` }}>
            <a
              href={`#${id}`}
              className={cn(
                "block text-sm hover:text-monokai-yellow transition-colors",
                activeId === id ? "text-monokai-green" : "text-monokai-text opacity-80",
              )}
            >
              {title}
            </a>
          </li>
        ))}
      </ul>
    </nav>
  )
}

