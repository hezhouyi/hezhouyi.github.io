# Public files at the repository root (CV PDF, profile.jpg, image/, code
# folders, QR codes in content/, …) are published at the same paths as in the
# repository. They are added as Jekyll static files, so `jekyll build` and
# `jekyll serve` both include them and keep them across rebuilds; only changed
# files are copied again. Jekyll's own pages win on a name clash.
require "set"

module RootFiles
  ROOT = File.expand_path("../..", __dir__)

  # Never published: tooling, sources and drafts.
  SKIP_TOP = %w[.git .github .gitignore .gitattributes site scripts draft README.md].freeze

  def self.skip?(rel)
    parts = rel.split("/")
    SKIP_TOP.include?(parts.first) ||
      parts.include?(".DS_Store") || parts.include?("__pycache__") ||
      (parts.first == "content" && rel.end_with?(".md"))
  end

  class Generator < Jekyll::Generator
    priority :lowest

    def generate(site)
      taken = (site.pages + site.static_files).map { |f| f.destination("").sub(%r{\A/}, "") }.to_set
      Dir.glob("**/*", File::FNM_DOTMATCH, base: ROOT).sort.each do |rel|
        next if rel.end_with?("/.", "/..") || rel == "." || RootFiles.skip?(rel)
        next unless File.file?(File.join(ROOT, rel))
        next if taken.include?(rel)

        site.static_files << Jekyll::StaticFile.new(site, ROOT, "/#{File.dirname(rel).sub(/\A\.\z/, "")}", File.basename(rel))
      end
    end
  end
end
