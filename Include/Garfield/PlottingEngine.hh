#ifndef G_PLOTTING_ENGINE_H
#define G_PLOTTING_ENGINE_H
#include <memory>

class TStyle;

namespace Garfield {

/// Plotting style.

class PlottingEngine {
 public:
  /// Default constructor.
  PlottingEngine();
  /// Destructor
  ~PlottingEngine();

  /// Use serif font.
  static void SetSerif();
  /// Use sans-serif font.
  static void SetSansSerif();

  /// Set the colour palette.
  static void SetPalette(const int ncol);

 private:
  static void SetFont(const int font);
  static std::unique_ptr<TStyle> m_style;
};

}  // namespace Garfield

#endif
