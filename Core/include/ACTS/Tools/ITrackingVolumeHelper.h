///////////////////////////////////////////////////////////////////
// ITrackingVolumeHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEHELPERH
#define ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEHELPERH 1

// Geometry module
#include "ACTS/Utilities/BinningType.h"
// Core module
#include <string>
#include <vector>
#include <memory>
#include "ACTS/Utilities/Definitions.h"

namespace Acts {

  class Layer;
  class TrackingVolume;
  class VolumeBounds;
  class Material;
  
  typedef std::shared_ptr<const Layer>               LayerPtr;
  typedef std::shared_ptr<const TrackingVolume>      TrackingVolumePtr;
  typedef std::shared_ptr<const VolumeBounds>        VolumeBoundsPtr;

  typedef std::vector< LayerPtr >                    LayerVector;
  typedef std::vector< TrackingVolumePtr >           TrackingVolumeVector;
  
  /** @class ITrackingVolumeHelper
  
     Interface class ITrackingVolumeHelper tools, it inherits from IAlgTool.
     The ITrackingVolumeHelper is a tool to pack a set of layers into a volume,
     or - to wrap several volumes into a container volume.
  
    TrackingVolumes only exist as std::shared_ptr

     @author Andreas.Salzburger@cern.ch
    */
  class ITrackingVolumeHelper
  {
    public:
    /**Virtual destructor*/
    virtual ~ITrackingVolumeHelper() = default;

    /** create a TrackingVolume* from a set of layers and (optional) parameters
          
	@param layers : vector of static layers confined by the TrackingVolume
	if no bounds or HepTransform is given, they define the size
	together with the volume enevlope parameters
	@param matprop : dense material properties for this TrackingVolume 
	@param volBounds : (optional) bounds of this TrackingVolume - ownership given
	@param transform : (optional) placement of this TrackingVolume - ownership given
	@param entryLayers : switch to build entry layers 
	@param volumeName  : volume name to be given
    */
    virtual TrackingVolumePtr createTrackingVolume(const LayerVector& layers,
						   const Material& matprop,
						   VolumeBoundsPtr volBounds,
						   std::shared_ptr<Transform3D> transform = nullptr,
						   const std::string& volumeName = "UndefinedVolume",
						   BinningType btype = arbitrary) const = 0;
                                                                                                            
    /** create a TrackingVolume* from a set of layers and (optional) parameters

	@param layers : vector of static layers confined by the TrackingVolume
	if no bounds or HepTransform is given, they define the size
	together with the volume enevlope parameters
	@param matprop : dense material properties for this TrackingVolume
	@param loc1Min, loc1Max, loc2Min, loc2Max : local position in space,
	this TrackingVolume is restricted to Translation only
	@param volumeName  : volume name to be given
    */
    virtual TrackingVolumePtr createTrackingVolume(const LayerVector& layers,
						   const Material& matprop,
						   double loc1Min, double loc1Max,
						   double loc2Min, double loc2Max,
						   const std::string& volumeName = "UndefinedVolume",
						   BinningType btype = arbitrary) const = 0;
                                                                                                            

    /** create a gap volume from dimensions and
       
	@param matprop : dense material properties for this TrackingVolume
	@param loc1Min, loc1Max, loc2Min, loc2Max : local position in space,
	this TrackingVolume is restricted to Translation only
	@param materialLayers : number of material layers (aequidistant binning)
	@param cylinder : type of layers 
	@param volumeName  : volume name to be given
                   
    */                                                      
    virtual TrackingVolumePtr createGapTrackingVolume(const Material& matprop,
						      double loc1Min, double loc1Max,
						      double loc2Min, double loc2Max,
						      unsigned int materialLayers,
						      bool cylinder = true,
						      const std::string& volumeName = "UndefinedVolume") const = 0;

    /** create a gap volume from dimensions and
       
	@param matprop : dense material properties for this TrackingVolume
	@param layerPositions : custom layer positions
	@param materialLayers : number of material layers (aequidistant binning)
	@param cylinder : type of layers 
	@param volumeName  : volume name to be given
                   
    */                                                      
    virtual TrackingVolumePtr createGapTrackingVolume(const Material& matprop,
						      double loc1Min, double loc1Max,
						      double loc2Min, double loc2Max,
						      const std::vector<double>& layerPositions,
						      bool cylinder = true,
						      const std::string& volumeName = "UndefinedVolume",
						      BinningType btype = arbitrary) const = 0;
                                                                                                                                                                  
    /** Create a one level higher TrackingVolue

	@param volumes : the volumes to be comnbined
	@param matprop : dense material properties for this TrackingVolume
	@param volumeName  : volume name to be given

    */
    virtual TrackingVolumePtr createContainerTrackingVolume(const TrackingVolumeVector& volumes) const = 0;
                                                                                                                                                         

  };

} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMECREATOR_H