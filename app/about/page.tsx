import Image from "next/image"
import { Github, X, Linkedin } from "lucide-react"

export default function About() {
  return (
    <div className="flex flex-col items-center">
      <div className="relative w-48 h-48 mb-8">
        <Image
          src="/resources/profile.jpeg"
          alt="Profile Picture"
          fill
          className="rounded-full object-cover"
        />
      </div>
      <h1 className="text-3xl font-bold mb-4 text-monokai-pink">About Me</h1>
      <p className="text-monokai-text max-w-2xl text-center mb-8">
        Hello! I'm a software developer currently working at Uber. I hold both a Bachelor's and Master's degree in Aerospace Engineering. My interests include programming, computational fluid dynamics and robotics. I started this blog to document new things I learn along the way.
        You can reach out to me at <a
          href="mailto:rohit@aerorohit.in"
          className="text-monokai-yellow hover:text-monokai-text transition-colors"
          aria-label="Email me at rohit@aerorohit.in"
        >
          rohit@aerorohit.in
        </a>.
      </p>
      <div className="flex space-x-6">
        <a
          href="https://x.com/aerorohit73"
          className="text-monokai-text hover:text-monokai-yellow transition-colors"
          aria-label="X"
        >
          <X size={24} />
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

