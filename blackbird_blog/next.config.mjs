/** @type {import('next').NextConfig} */
const nextConfig = {

    output: 'export', // Enable static export
    trailingSlash: true, // Recommended for GitHub Pages routing
    // If using a project repository (<username>.github.io/<repo-name>)
    basePath: process.env.NODE_ENV === 'production' ? '/aerorohit.github.io' : '',
    assetPrefix: process.env.NODE_ENV === 'production' ? '/aerorohit.github.io' : '',
    images: {
        unoptimized: true,
    },
};

export default nextConfig;
