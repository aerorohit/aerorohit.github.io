import Image from "next/image"
import { Github, X as XIcon, Linkedin } from "lucide-react"

export default function About() {
  return (
    <div className="flex flex-col items-center">
      <div className="relative w-48 h-48 mb-8">
        <Image
          src="/resources/rohit_blog_about.jpeg"
          alt="Profile Picture"
          fill
          className="rounded-full object-cover"
        />
      </div>
      <h1 className="text-3xl font-bold mb-4 text-monokai-pink">About Me</h1>
      <p className="text-monokai-text max-w-2xl text-center mb-8">
      Jotting down thoughts as they come. I am a software developer at Uber with interests in programming and computational fluid dynamics. I am also exploring robotics and AI.
      </p>
      <div className="flex space-x-6">
        <a
          href="https://x.com/aerorohit73"
          className="text-monokai-text hover:text-monokai-yellow transition-colors"
          aria-label="X"
        >
          <XIcon size={24} />
        </a>
        <a
          href="https://github.com/aerorohit"
          className="text-monokai-text hover:text-monokai-yellow transition-colors"
          aria-label="GitHub"
        >
          <Github size={24} />
        </a>
        <a
          href="https://linkedin.com/in/aerorohit"
          className="text-monokai-text hover:text-monokai-yellow transition-colors"
          aria-label="LinkedIn"
        >
          <Linkedin size={24} />
        </a>
      </div>
    </div>
  )
}
