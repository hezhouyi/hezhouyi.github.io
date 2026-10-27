# The editable content of the site lives in content/ at the repository root,
# one Markdown file per page or list: settings (YAML front matter) on top,
# free text below. This plugin makes each file available to the templates as
#   site.data.content.<file name>          e.g. site.data.content.about.status
# with the text below the front matter as `.body`.
# scripts/check_content.py checks these files (with clear messages) before Jekyll runs.
#
# Links and file locations in content/ follow one rule: a full web address
# (https://…, mailto:…), or a file or page of this site written from the top
# folder of the repository, starting with "/" and spelled as on GitHub (spaces
# allowed), e.g. "/profile.jpg". Templates pass them through `link`.
require "date"
require "kramdown"
require "open3"
require "yaml"

module ContentFolder
  DIR = File.expand_path("../../content", __dir__)

  def self.read(path)
    text = File.read(path, encoding: "utf-8")
    if (m = text.match(/\A---\s*\n(.*?)\n---[ \t]*\n?(.*)\z/m))
      data = YAML.safe_load(m[1], permitted_classes: [Date, Time]) || {}
      [data, m[2]]
    else
      [{}, text]
    end
  rescue Psych::SyntaxError => e
    raise Jekyll::Errors::FatalException, "#{path}: #{e.message}"
  end

  # Empty settings ("", only spaces, empty lists) are dropped, so templates can
  # simply write {% if a.greeting %} and a left-out setting leaves no trace.
  def self.clean(value)
    case value
    when String then value.strip.empty? ? nil : value
    when Array then value.map { |v| clean(v) }.compact.then { |a| a.empty? ? nil : a }
    when Hash then value.transform_values { |v| clean(v) }.compact.then { |h| h.empty? ? nil : h }
    else value
    end
  end

  # "/Zhouyi He.pdf" -> "/Zhouyi%20He.pdf". Already encoded parts stay as they are,
  # so this can be applied more than once; what follows ? or # is left alone.
  UNSAFE = %r{[^A-Za-z0-9\-._~!$&'()*+,;=:@/]}

  def self.encode_path(href)
    path, rest = href.split(/(?=[?#])/, 2)
    path = path.b.gsub(/%\h\h/n) { |m| m[1, 2].hex.chr }
    path.gsub(UNSAFE) { |c| format("%%%02X", c.ord) }.force_encoding("UTF-8") + rest.to_s
  end

  class Generator < Jekyll::Generator
    priority :highest

    def generate(site)
      site.data["content"] = {}
      Dir[File.join(DIR, "*.md")].sort.each do |path|
        name = File.basename(path, ".md")
        data, body = ContentFolder.read(path)
        site.data["content"][name] = ContentFolder.clean(data.merge("body" => body)) || {}
      end
      # "Last updated" on the Resume page: when the CV PDF was last committed. Looked up
      # on every build, so `jekyll serve` shows the same date as the published site.
      site.data["cv_updated"] = ContentFolder.last_commit_date(site.data.dig("content", "cv", "pdf"))
    end
  end

  # The date (YYYY-MM-DD) of the last commit that changed a file given as in
  # content/ ("/Zhouyi He.pdf"); nil when git or the file's history is not available.
  def self.last_commit_date(location)
    return nil unless location
    path = location.to_s.b.gsub(/%\h\h/n) { |m| m[1, 2].hex.chr }.force_encoding("UTF-8").delete_prefix("/")
    out, status = Open3.capture2("git", "log", "-1", "--format=%cs", "--", path, chdir: File.dirname(DIR))
    date = out.strip
    status.success? && !date.empty? ? date : nil
  rescue SystemCallError
    nil
  end

  # Liquid filters for short texts written in Markdown inside the front matter.
  # Quotes are left exactly as typed.
  module Filters
    OPTIONS = { input: "kramdown", smart_quotes: %w[apos apos quot quot], auto_ids: false }.freeze

    # A paper's summary (publications.md) shows at most this many words; a longer
    # one is cut there and ends in "…". Also stated in that file's comment.
    SUMMARY_MAX_WORDS = 50

    # Like inline_md, cut after SUMMARY_MAX_WORDS words. Only the visible text is
    # counted, and links or italics open at the cut are closed properly.
    def summary_md(text)
      words = 0
      open_tags = []
      out = +""
      inline_md(text).scan(%r{<[^>]+>|[^<]+}) do |token|
        if token.start_with?("<")
          if token.start_with?("</") then open_tags.pop
          elsif !token.end_with?("/>") && !token.start_with?("<br") then open_tags.push(token[/\A<(\w+)/, 1])
          end
          out << token
          next
        end
        token.split(/(\s+)/).each do |piece|
          words += 1 unless piece.strip.empty?
          if words > SUMMARY_MAX_WORDS
            return out.rstrip + "…" + open_tags.reverse.map { |t| "</#{t}>" }.join
          end
          out << piece
        end
      end
      out
    end

    # A one-line text: links, bold, italic — without a surrounding paragraph.
    def inline_md(text)
      html = Kramdown::Document.new(text.to_s, OPTIONS).to_html.strip
      html.sub(%r{\A<p>(.*)</p>\z}m, '\1').gsub("\n", "<br>\n")
    end

    # Paragraphs, where lines written under each other stay on separate lines.
    def lines_md(text)
      Kramdown::Document.new(text.to_s, OPTIONS).to_html.gsub(%r{<p>(.*?)</p>}m) do
        "<p>#{Regexp.last_match(1).strip.gsub("\n", "<br>\n")}</p>"
      end
    end

    # A link or file location from content/ as a URL; follow with relative_url.
    def link(value)
      value = value.to_s.strip
      value.start_with?("/") ? ContentFolder.encode_path(value) : value
    end
  end
end

Liquid::Template.register_filter(ContentFolder::Filters)

# Links to other websites open in a new tab, wherever they come from; links to
# this site written in Markdown text ([text](/Zhouyi He.pdf)) get their spaces
# encoded, and the site's path prefix when it is not at the root (a preview).
Jekyll::Hooks.register :pages, :post_render do |page|
  next unless page.output_ext == ".html"

  base = page.site.baseurl.to_s.chomp("/")
  page.output = page.output.gsub(%r{(<a href="|<img src=")(/[^"]*)"}) do
    path = Regexp.last_match(2)
    path = base + path unless base.empty? || path == base || path.start_with?("#{base}/")
    %(#{Regexp.last_match(1)}#{ContentFolder.encode_path(path)}")
  end
  page.output = page.output.gsub(/<a href="(https?:[^"]*)"(?![^>]*\btarget=)/) do
    %(<a href="#{Regexp.last_match(1)}" target="_blank" rel="noopener")
  end
end
