# Structured data for search engines (schema.org Person), built from
# content/seo.md ("person:") and content/site.md (the name). Settings that are
# left out are simply not included.
require "json"

module SeoFilters
  def person_jsonld(content)
    site = @context.registers[:site]
    person = content.dig("seo", "person") || {}
    orgs = lambda do |list, type|
      list&.map { |o| { "@type" => o["type"] || type, "name" => o["name"], "url" => o["url"] }.compact }
    end
    data = {
      "@context" => "https://schema.org",
      "@type" => "Person",
      "name" => content.dig("site", "name"),
      "alternateName" => content.dig("site", "name_cn"),
      "url" => "#{site.config["url"]}#{site.baseurl}/",
      "image" => (site.data.dig("generated", "portrait") && "#{site.config["url"]}#{site.baseurl}/assets/img/portrait.jpg"),
      "jobTitle" => person["job_title"],
      "description" => person["description"],
      "affiliation" => orgs.call(person["works_at"], "Organization"),
      "alumniOf" => orgs.call(person["studied_at"], "CollegeOrUniversity"),
      "knowsAbout" => person["knows_about"],
      "sameAs" => person["profiles"],
    }.compact
    # "</" could end the <script> element early
    JSON.pretty_generate(data).gsub("</", "<\\/")
  end
end

Liquid::Template.register_filter(SeoFilters)
