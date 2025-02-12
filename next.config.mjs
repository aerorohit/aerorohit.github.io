/** @type {import('next').NextConfig} */
const nextConfig = {

    output: 'export', // Enable static export
    trailingSlash: true, // Recommended for GitHub Pages routing
    // If using a project repository (<username>.github.io/<repo-name>)
    // basePath: 'aerorohit.github.io',
    assetPrefix: '',
};

export default nextConfig;
