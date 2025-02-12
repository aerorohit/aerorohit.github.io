import Image from "next/image"
import { Github, Twitter, Linkedin } from "lucide-react"

export default function About() {
  return (
    <div className="flex flex-col items-center">
      <div className="relative w-48 h-48 mb-8">
        <Image
          src="/resources/profile.png"
          alt="Profile Picture"
          fill
          className="rounded-full object-cover"
        />
      </div>
      <h1 className="text-3xl font-bold mb-4 text-monokai-pink">About Me</h1>
      <p className="text-monokai-text max-w-2xl text-center mb-8">
        Hello! I'm a passionate developer and writer. I love exploring new technologies and sharing my knowledge through
        this blog. When I'm not coding, you can find me reading sci-fi novels or hiking in the great outdoors.
      </p>
      <div className="flex space-x-6">
        <a
          href="https://twitter.com/yourusername"
          className="text-monokai-text hover:text-monokai-yellow transition-colors"
          aria-label="Twitter"
        >
          <Twitter size={24} />
        </a>
        <a
          href="https://github.com/yourusername"
          className="text-monokai-text hover:text-monokai-yellow transition-colors"
          aria-label="GitHub"
        >
          <Github size={24} />
        </a>
        <a
          href="https://linkedin.com/in/yourusername"
          className="text-monokai-text hover:text-monokai-yellow transition-colors"
          aria-label="LinkedIn"
        >
          <Linkedin size={24} />
        </a>
      </div>
    </div>
  )
}

