// Photography page: open a photo in the lightbox, step through with the arrows or keyboard.
(function () {
  const items = Array.from(document.querySelectorAll(".gallery button"));
  const lb = document.getElementById("lightbox");
  if (!items.length || !lb) return;

  const lbImg = lb.querySelector("img");
  const lbText = lb.querySelector(".caption span");
  const lbLink = lb.querySelector(".caption a");
  let current = 0;

  function show(i) {
    current = (i + items.length) % items.length;
    const d = items[current].dataset;
    lbImg.src = d.large;
    lbImg.alt = d.caption;
    lbText.textContent = d.caption;
    lbLink.href = d.original;
  }
  function openLb(i) { show(i); lb.classList.add("open"); document.body.style.overflow = "hidden"; }
  function closeLb() { lb.classList.remove("open"); document.body.style.overflow = ""; items[current].focus(); }

  items.forEach((b, i) => b.addEventListener("click", () => openLb(i)));
  lb.querySelector(".lb-close").onclick = closeLb;
  lb.querySelector(".lb-prev").onclick = () => show(current - 1);
  lb.querySelector(".lb-next").onclick = () => show(current + 1);
  lb.addEventListener("click", e => { if (e.target === lb) closeLb(); });
  document.addEventListener("keydown", e => {
    if (!lb.classList.contains("open")) return;
    if (e.key === "Escape") closeLb();
    if (e.key === "ArrowLeft") show(current - 1);
    if (e.key === "ArrowRight") show(current + 1);
  });
})();
