# The build system for Adobe's "XMP Toolkit SDK" is just absolutely bonkers.
# Here we apply our own build system to a submodule checkout. It just generates
# a static library that our main executable links to.

cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(EXPAT REQUIRED)
find_package(ZLIB REQUIRED)

include_directories(
    .
    ${CMAKE_CURRENT_SOURCE_DIR}/XMP-Toolkit-SDK
    ${CMAKE_CURRENT_SOURCE_DIR}/XMP-Toolkit-SDK/public/include
    ${CMAKE_CURRENT_SOURCE_DIR}/XMP-Toolkit-SDK/XMPFilesPlugins/api/source
    ${EXPAT_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
)

add_compile_definitions(
    UNIX_ENV=1
    XMP_UNIXBuild=1
    XMP_COMPONENT_INT_NAMESPACE=AdobeXMPCore_Int
)

add_custom_target(
    expat_header_dir
    ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory third-party/expat/lib
)

add_custom_target(
    expat_header_link
    ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${EXPAT_INCLUDE_DIRS}/expat.h third-party/expat/lib/expat.h
)

add_dependencies(expat_header_link expat_header_dir)

add_custom_target(
    zlib_header_dir
    ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory third-party/zlib
)

add_custom_target(
    zlib_header_link
    ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${ZLIB_INCLUDE_DIRS}/zlib.h third-party/zlib/zlib.h
)

add_dependencies(zlib_header_link zlib_header_dir)

add_library(
    adobexmp
    STATIC
    XMP-Toolkit-SDK/source/EndianUtils.hpp
    XMP-Toolkit-SDK/source/ExpatAdapter.hpp
    XMP-Toolkit-SDK/source/Host_IO.hpp
    XMP-Toolkit-SDK/source/Host_IO-POSIX.cpp
    XMP-Toolkit-SDK/source/IOUtils.cpp
    XMP-Toolkit-SDK/source/IOUtils.hpp
    XMP-Toolkit-SDK/source/PerfUtils.cpp
    XMP-Toolkit-SDK/source/PerfUtils.hpp
    XMP-Toolkit-SDK/source/SafeStringAPIs.cpp
    XMP-Toolkit-SDK/source/UnicodeConversions.cpp
    XMP-Toolkit-SDK/source/UnicodeConversions.hpp
    XMP-Toolkit-SDK/source/XIO.cpp
    XMP-Toolkit-SDK/source/XIO.hpp
    XMP-Toolkit-SDK/source/XML_Node.cpp
    XMP-Toolkit-SDK/source/XMLParserAdapter.hpp
    XMP-Toolkit-SDK/source/XMPFiles_IO.cpp
    XMP-Toolkit-SDK/source/XMPFiles_IO.hpp
    XMP-Toolkit-SDK/source/XMP_LibUtils.cpp
    XMP-Toolkit-SDK/source/XMP_LibUtils.hpp
    XMP-Toolkit-SDK/source/XMP_ProgressTracker.cpp
    XMP-Toolkit-SDK/source/XMP_ProgressTracker.hpp
    XMP-Toolkit-SDK/third-party/zuid/interfaces/MD5.cpp
    XMP-Toolkit-SDK/XMPCommon/source/AutoSharedLock.cpp
    XMP-Toolkit-SDK/XMPCommon/source/ConfigurableImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/ConfigurationManagerImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/DefaultErrorNotifierImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/ErrorImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/ErrorNotifierWrapperImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/IConfigurable_I.cpp
    XMP-Toolkit-SDK/XMPCommon/source/IConfigurationManager_I.cpp
    XMP-Toolkit-SDK/XMPCommon/source/IError_I.cpp
    XMP-Toolkit-SDK/XMPCommon/source/IErrorNotifier_I.cpp
    XMP-Toolkit-SDK/XMPCommon/source/IMemoryAllocator_I.cpp
    XMP-Toolkit-SDK/XMPCommon/source/IUTF8StringComparator.cpp
    XMP-Toolkit-SDK/XMPCommon/source/IUTF8String_I.cpp
    XMP-Toolkit-SDK/XMPCommon/source/MemoryAllocatorWrapperImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/MemoryManagedObject.cpp
    XMP-Toolkit-SDK/XMPCommon/source/SharedMutexImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/SharedObjectImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/ThreadSafeImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/UTF8StringImpl.cpp
    XMP-Toolkit-SDK/XMPCommon/source/XMPCommon.cpp
    XMP-Toolkit-SDK/XMPCore/source/ArrayNodeImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/ClientDOMParserWrapperImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/ClientDOMSerializerWrapperImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/CompositeNodeImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/CoreConfigurationManagerImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/CoreObjectFactoryImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/DOMImplementationRegistryImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/DOMParserImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/DOMSerializerImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/ExpatAdapter.cpp
    XMP-Toolkit-SDK/XMPCore/source/IArrayNode_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/ICompositeNode_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/ICoreConfigurationManager_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/ICoreObjectFactory_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IDOMImplementationRegistry_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IDOMParser_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IDOMSerializer_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IMetadataConverterUtils_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IMetadata_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/INameSpacePrefixMap_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/INode_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/INodeIterator_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IPath_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IPathSegment_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/ISimpleNode_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/IStructureNode_I.cpp
    XMP-Toolkit-SDK/XMPCore/source/MetadataConverterUtilsImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/MetadataImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/NameSpacePrefixMapImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/NodeImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/ParseRDF.cpp
    XMP-Toolkit-SDK/XMPCore/source/PathImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/PathSegmentImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/RDFDOMParserImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/RDFDOMSerializerImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/SimpleNodeImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/StructureNodeImpl.cpp
    XMP-Toolkit-SDK/XMPCore/source/WXMPIterator.cpp
    XMP-Toolkit-SDK/XMPCore/source/WXMPMeta.cpp
    XMP-Toolkit-SDK/XMPCore/source/WXMPUtils.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPCore_Impl.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPCore_Impl.hpp
    XMP-Toolkit-SDK/XMPCore/source/XMPIterator2.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPIterator2.hpp
    XMP-Toolkit-SDK/XMPCore/source/XMPIterator.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPIterator.hpp
    XMP-Toolkit-SDK/XMPCore/source/XMPMeta2-GetSet.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPMeta2.hpp
    XMP-Toolkit-SDK/XMPCore/source/XMPMeta.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPMeta-GetSet.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPMeta.hpp
    XMP-Toolkit-SDK/XMPCore/source/XMPMeta-Parse.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPMeta-Serialize.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPUtils2.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPUtils.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPUtils-FileInfo.cpp
    XMP-Toolkit-SDK/XMPCore/source/XMPUtils.hpp
    XMP-Toolkit-SDK/XMPFilesPlugins/api/source/HostAPIAccess.cpp
    XMP-Toolkit-SDK/XMPFilesPlugins/api/source/PluginAPIImpl.cpp
    XMP-Toolkit-SDK/XMPFilesPlugins/api/source/PluginBase.cpp
    XMP-Toolkit-SDK/XMPFilesPlugins/api/source/PluginRegistry.cpp
    XMP-Toolkit-SDK/XMPFilesPlugins/api/source/PluginUtils.cpp
    XMP-Toolkit-SDK/XMPFilesPlugins/PluginTemplate/source/Template_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/AIFF_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/AIFF_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/ASF_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/ASF_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/AVCHD_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/AVCHD_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/Basic_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/Basic_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/FLV_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/FLV_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/GIF_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/GIF_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/InDesign_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/InDesign_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/JPEG_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/JPEG_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/MP3_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/MP3_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/MPEG2_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/MPEG2_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/MPEG4_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/MPEG4_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/P2_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/P2_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/PNG_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/PNG_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/PostScript_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/PostScript_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/PSD_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/PSD_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/RIFF_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/RIFF_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/Scanner_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/Scanner_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/SonyHDV_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/SonyHDV_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/SVG_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/SVG_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/SWF_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/SWF_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/TIFF_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/TIFF_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/Trivial_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/Trivial_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/UCF_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/UCF_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/WAVE_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/WAVE_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAMEX_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAMEX_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAMFAM_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAMFAM_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAM_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAM_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAMSAM_Handler.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FileHandlers/XDCAMSAM_Handler.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/AIFF/AIFFBehavior.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/AIFF/AIFFMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/AIFF/AIFFReconcile.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ASF_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ASF_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ID3_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ID3_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/IFF/ChunkController.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/IFF/Chunk.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/IFF/ChunkPath.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/IFF/IChunkBehavior.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/IPTC_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/IPTC_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ISOBaseMedia_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ISOBaseMedia_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/META_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/META_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/MOOV_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/MOOV_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/P2_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/P2_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PackageFormat_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PackageFormat_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PNG_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PNG_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PostScript_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PostScript_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PSIR_FileWriter.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PSIR_MemoryReader.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/PSIR_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/QuickTime_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/QuickTime_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/Reconcile_Impl.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/Reconcile_Impl.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ReconcileIPTC.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ReconcileLegacy.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ReconcileLegacy.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/ReconcileTIFF.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/RIFF.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/RIFF.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/RIFF_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/RIFF_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/SVG_Adapter.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/SVG_Adapter.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/SWF_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/SWF_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/TIFF_FileWriter.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/TIFF_MemoryReader.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/TIFF_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/TIFF_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/TimeConversionUtils.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/TimeConversionUtils.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/BEXTMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/CartMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/Cr8rMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/DISPMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/INFOMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/iXMLMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/PrmLMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/WAVEBehavior.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/WAVE/WAVEReconcile.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/XDCAM_Support.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/XDCAM_Support.hpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/XMPScanner.cpp
    XMP-Toolkit-SDK/XMPFiles/source/FormatSupport/XMPScanner.hpp
    XMP-Toolkit-SDK/XMPFiles/source/HandlerRegistry.cpp
    XMP-Toolkit-SDK/XMPFiles/source/NativeMetadataSupport/IMetadata.cpp
    XMP-Toolkit-SDK/XMPFiles/source/NativeMetadataSupport/IReconcile.cpp
    XMP-Toolkit-SDK/XMPFiles/source/NativeMetadataSupport/MetadataSet.cpp
    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/FileHandlerInstance.cpp
    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/HostAPIImpl.cpp
    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/Module.cpp
    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/OS_Utils_Linux.cpp
    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/PluginManager.cpp
    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/XMPAtoms.cpp
    XMP-Toolkit-SDK/XMPFiles/source/WXMPFiles.cpp
    XMP-Toolkit-SDK/XMPFiles/source/XMPFiles.cpp
    XMP-Toolkit-SDK/XMPFiles/source/XMPFiles.hpp
    XMP-Toolkit-SDK/XMPFiles/source/XMPFiles_Impl.cpp
    XMP-Toolkit-SDK/XMPFiles/source/XMPFiles_Impl.hpp
)

add_dependencies(adobexmp expat_header_link zlib_header_link)

#    XMP-Toolkit-SDK/source/Host_IO-Win.cpp
#    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/OS_Utils_Android.cpp
#    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/OS_Utils_Mac.cpp
#    XMP-Toolkit-SDK/XMPFiles/source/PluginHandler/OS_Utils_WIN.cpp
