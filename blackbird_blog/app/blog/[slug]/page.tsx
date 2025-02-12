import fs from 'fs';
import path from 'path';
import { getBlogPost} from "@/lib/mdx"
import { MDXRemote } from "next-mdx-remote/rsc"
import rehypeHighlight from "rehype-highlight"
import rehypeKatex from "rehype-katex"
import remarkMath from "remark-math"
import { formatDate } from "@/lib/utils"
// import TableOfContents from "@/components/TableOfContents"
import 'katex/dist/katex.min.css'
import 'highlight.js/styles/monokai-sublime.css'

const postsDirectory = path.join(process.cwd(), 'posts'); 
export async function generateStaticParams() {
  const filenames = fs.readdirSync(postsDirectory);

  return filenames.map((filename) => ({
    slug: filename.replace(/\.mdx$/, ''),
  }));
}

export default async function BlogPost({ params }: { params: { slug: string } }) {
  const post = await getBlogPost(params.slug)
  // const toc = getTableOfContents(post.content)

  return (
    <div className="flex gap-8">
      <article className="prose prose-invert mx-auto flex-grow max-w-3xl">
        <header className="mb-8 pb-8 border-b border-monokai-lighter">
          <h1 className="text-4xl font-bold mb-4 text-monokai-orange">{post.title}</h1>
          <time className="text-monokai-text opacity-60">{formatDate(post.date)}</time>
        </header>
        <MDXRemote
          source={post.content}
          options={{
            mdxOptions: {
              remarkPlugins: [remarkMath],
              rehypePlugins: [rehypeHighlight, rehypeKatex],
            },
          }}
        />
      </article>
      {/* <aside className="hidden lg:block w-64 flex-shrink-0">
        <TableOfContents items={toc} />
      </aside> */}
    </div>
  )
}

