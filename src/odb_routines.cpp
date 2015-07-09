// System includes
#include <iostream>
using namespace std;
// Begin local includes
#include <odb_API.h>
#include <odb_MaterialTypes.h>
#include <odb_SectionTypes.h>
// To do exceptions
extern "C"
{
//****************************************************************************80
// Model data
//****************************************************************************80
  odb_Odb *model_data(  char *odb_name,
                        char *fe_type_name,
                        int dom,
                        int nnod,
                        float *coo,
                        int nel,
                        int nnel,
                        int *conn)
  {
    int i, j;
    odb_Odb *odbPtr;
    // Create an ODB (which also creates the rootAssembly).
    odb_String odbName(odb_name);
    odb_Odb& odb = Odb("simpleModel",
                       "ODB created with C++ ODB API",
                       "example illustrating C++ ODB API",
                       odbName);
      // Section and part data
    odb_SectionApi sectionApi;
    odb.extendApi(odb_Enum::odb_SECTION,sectionApi);
    odb_String sectionName("Homogeneous Solid Section");
    odb_String materialName("Elastic material");
    odb_HomogeneousSolidSection& section_1 =
    sectionApi.HomogeneousSolidSection(sectionName,materialName);
    odb_Part part1;
    if ( dom == 2 ){
      part1 = odb.Part("part-1",odb_Enum::TWO_D_PLANAR,
                       odb_Enum::DEFORMABLE_BODY);
    }
    else if ( dom == 3){
      part1 = odb.Part("part-1",odb_Enum::THREE_D,
                       odb_Enum::DEFORMABLE_BODY);
    }
    // Node data
    odb_SequenceInt nodeLabels;
    for(i=1; i<nnod+1; i++)
      nodeLabels.append(i);
    odb_SequenceSequenceFloat nodeCoor;
    for(i=0; i<nnod; i++) {
      odb_SequenceFloat loc;
      for(j=0; j<dom; j++) {
        loc.append(coo[nnod*j+i]);
      }
      nodeCoor.append(loc);
    }
    odb_String nodeSetName("nset-1");
    part1.addNodes(nodeLabels,nodeCoor,nodeSetName);
    // Element data
    odb_SequenceInt elLabels;
    for(i=1; i<nel+1; i++)
      elLabels.append(i);
    odb_SequenceSequenceInt elConnect;
    for(i=0; i<nel; i++) {
      odb_SequenceInt loc;
      for(j=0; j<nnel; j++) {
        loc.append(conn[nel*j+i]);
      }
      elConnect.append(loc);
    }
    odb_String elType(fe_type_name);
    odb_String elsetName("eset-1");
    part1.addElements(elLabels,elConnect,elType,elsetName);
    // Instance the part
    odb_String partInstanceName("part-1-1");
    odb_Instance& instance1 = odb.rootAssembly().Instance(partInstanceName, part1);
    // Create instance level sets for section assignment
    elsetName = "Material 1";
    odb_Set& elset_1 = instance1.ElementSet(elsetName,elLabels);
    // Section assignment on instance
    instance1.assignSection(elset_1,section_1);
    // Save and point to odb
    odb.save();
    odbPtr = &odb;
    return odbPtr;
  }
//****************************************************************************80
// Close odb
//****************************************************************************80
  void close_odb(odb_Odb &odb)
  {
    odb.close();
  }
//****************************************************************************80
// Create step
//****************************************************************************80
  void step(odb_Odb &odb,
            char *step_name,
            char *step_desc )
  {
    odb_String stepName(step_name);
    odb_String stepDescription(step_desc);
    odb_Step &step1 = odb.Step(stepName,
                               stepDescription,
                               odb_Enum::TIME,1.0);
    // Save
    odb.save();
  }

  // Create frame
  void frame(odb_Odb &odb,
             char* step_name,
             char* frame_desc,
             int fnum,
             float analy_time)
  {
    // Open step
    odb_String stepName(step_name);
    odb_Step& step1 = odb.steps()[stepName];
    // Create frame
    odb_String frameDescription(frame_desc);
    odb_Frame frame1 = step1.Frame(fnum,
                                   analy_time,
                                   frameDescription);
    // Save
    odb.save();
  }

  // Wite scalar field
  void scalar_field(odb_Odb &odb,
                     char *step_name,
                     int fnum,
                     char *field_name,
                     char *field_desc,
                     int nnod,
                     float* field_vals,
                     bool df)
  {
    int i;
    // Open step and frame
    odb_String stepName(step_name);
    odb_Step& step1 = odb.steps()[stepName];
    odb_Frame& frame1 = odb.steps()[stepName].frames(fnum-1);
    // Open instance
    odb_Instance& instance1 = odb.rootAssembly().instances()["part-1-1"];
    // Create field
    odb_String fieldName(field_name);
    odb_String fieldDescription(field_desc);
    odb_FieldOutput& scaField =
    frame1.FieldOutput(fieldName,
                       fieldDescription,
                       odb_Enum::SCALAR);
    // Get node sequence (bedarija, morm pogledat kuko se dobi vn stevilo
    // nodov v pythonu )
    odb_SequenceInt nodeLabels;
    for(i=1; i<nnod+1; i++){
      nodeLabels.append(i);
    }
    // Field data
    odb_SequenceSequenceFloat fieldData;
    for(i=0; i<nnod; i++){
      odb_SequenceFloat loc;
      loc.append(field_vals[i]);
      fieldData.append(loc);
    }
    // Add data
    scaField.addData(odb_Enum::NODAL,
                     instance1,
                     nodeLabels,
                     fieldData);
    // Make this the default field for visualization.
    if(df){
      step1.setDefaultField(scaField);
    }
    // Save
    odb.save();
  }

}
