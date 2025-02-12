import Link from "next/link"
import { getBlogPosts } from "@/lib/mdx"
import { formatDate } from "@/lib/utils"

export default async function Home() {
  const posts = await getBlogPosts()

  return (
    <div className="flex flex-col md:flex-row">
      <div className="md:w-1/3 mb-8 md:mb-0 md:pr-8">
        <h1 className="text-3xl font-bold mb-4 text-monokai-orange">Rohit Tembhare's website</h1>
        <p className="text-monokai-text">
        Sharing my thoughts as they come. I’m a software developer at Uber, passionate about programming and computational fluid dynamics, with an ongoing exploration of robotics and AI.
        </p>
      </div>
      <div className="md:border-l md:border-monokai-lighter md:pl-8 md:w-2/3">
        <h2 className="text-2xl font-bold mb-4 text-monokai-orange">Posts</h2>
        <ul className="space-y-4">
          {posts.map((post) => (
            <li key={post.slug}>
              <Link
                href={`/blog/${post.slug}`}
                className="block hover:bg-monokai-lighter p-4 rounded transition-colors"
              >
                <h3 className="text-xl font-semibold mb-2 text-monokai-green">{post.title}</h3>
                <time className="text-sm text-monokai-text opacity-60 mb-2 block">{formatDate(post.date)}</time>
                <p className="text-monokai-text opacity-80">{post.excerpt}</p>
              </Link>
            </li>
          ))}
        </ul>
      </div>
    </div>
  )
}

