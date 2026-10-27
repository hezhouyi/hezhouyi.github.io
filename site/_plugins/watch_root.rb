# `jekyll serve` only watches this folder (site/). This also watches content/
# and the public files at the repository root (see root_files.rb), so editing
# content/about.md or replacing the CV PDF rebuilds the preview by itself.
# Photos and the About-page photo are resized by Python, not here: after
# changing those, run scripts/build_gallery.py or scripts/build_portrait.py.
module WatchRoot
  # Folders whose changes never affect the site (and site/ is watched already)
  IGNORE = [%r{\A(site|scripts|draft|\.git|\.github)/}, /__pycache__/, /\.DS_Store\z/].freeze

  Jekyll::Hooks.register :site, :post_write do |site|
    next unless site.config["watch"] && !@listener

    require "listen"
    @listener = Listen.to(RootFiles::ROOT, ignore: IGNORE) do |modified, added, removed|
      changed = (modified + added + removed).map { |p| p.delete_prefix("#{RootFiles::ROOT}/") }
      Jekyll.logger.info "Regenerating:", "#{changed.first}#{" and #{changed.size - 1} more" if changed.size > 1}"
      begin
        site.process
        Jekyll.logger.info "", "done."
      rescue StandardError => e
        Jekyll.logger.error "Error:", e.message
      end
    end
    @listener.start
  end
end
