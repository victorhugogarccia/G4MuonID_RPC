#include "DetectorConstruction.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include <cmath>
#include <math.h>

#include "G4VSensitiveDetector.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

DetectorConstruction::DetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume *DetectorConstruction::Construct()
{


////////////////////////////////////////////////////////////
//Intro de materiales


	G4NistManager *nist = G4NistManager::Instance();
	
	//Mezcla para freon
	 G4Element* elC = nist->FindOrBuildElement("C");
    G4Element* elCl = nist->FindOrBuildElement("Cl");
    G4Element* elF = nist->FindOrBuildElement("F");
    G4Element* elH = nist->FindOrBuildElement("H");
    G4Element* elO = nist->FindOrBuildElement("O");
    // Definir las fracciones de masa para cada componente
    G4double C_fraction = 1;   // Número de átomos de carbono en el Freón-12
    G4double Cl_fraction = 2;  // Número de átomos de cloro en el Freón-12
    G4double F_fraction = 2;   // Número de átomos de flúor en el Freón-12

    // Calcular la fracción molar normalizada para cada componente
    G4double total_fraction = C_fraction + Cl_fraction + F_fraction;
    G4double C_molar_fraction = C_fraction / total_fraction;
    G4double Cl_molar_fraction = Cl_fraction / total_fraction;
    G4double F_molar_fraction = F_fraction / total_fraction;

    // Definir las densidades para cada componente
    G4double C_density = 2.26 * g/cm3;  // Densidad del carbono en el Freón-12
    G4double Cl_density = 1.56 * g/cm3; // Densidad del cloro en el Freón-12
    G4double F_density = 1.55 * g/cm3;  // Densidad del flúor en el Freón-12

    // Calcular la densidad ponderada de la mezcla de Freón
    G4double mixture_density = (C_density * C_fraction + Cl_density * Cl_fraction + F_density * F_fraction) / total_fraction;

    // Crear el material de la mezcla de Freón
    G4Material* freonMixture = new G4Material("FreonMixture", mixture_density, 3);
    freonMixture->AddElement(elC, C_molar_fraction);
    freonMixture->AddElement(elCl, Cl_molar_fraction);
    freonMixture->AddElement(elF, F_molar_fraction);

	//fin
	
	
	G4Material *steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");//PARA CREAR EL ACERO
	G4Material *glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");//Para el vidrio
	G4Material *grafito = nist->FindOrBuildMaterial("G4_GRAPHITE");

	G4Material *cobre = nist->FindOrBuildMaterial("G4_Cu");

// PMMA C5H8O2 ( Acrylic )
//  density = 1.19*g/cm3;
//  G4Material* Acrylic = new G4Material(name="Acrylic", density, numel=3);
//  Acrylic->AddElement(elC, 5);
//  Acrylic->AddElement(elH, 8);
//  Acrylic->AddElement(elO, 2);
	
//Para crear el material acrilico 
G4double density = 1.19*g/cm3;
	G4String materialName = "acrylic";
	G4Material* acrylic = new G4Material("acrylic", density, 3);
	acrylic->AddElement(elC, 5);
	acrylic->AddElement(elH, 8);
	acrylic->AddElement(elO, 2);
	
	//G4Material *SiO2 = new G4Material("SiO2", 2.201*g/cm3, 2);
	//SiO2->AddElement(nist->FindOrBuildElement("Si"), 1);
	//SiO2->AddElement(nist->FindOrBuildElement("O"), 2);
	
	//G4Material *H2O = new G4Material("H2O", 1.000*g/cm3, 2);
	//H2O->AddElement(nist->FindOrBuildElement("H"), 2);
	//H2O->AddElement(nist->FindOrBuildElement("O"), 1);
	
	//G4Element *C = nist->FindOrBuildElement("C");
	
	//G4Material *Aerogel = new G4Material("Aerogel", 0.200*g/cm3, 3);
	//Aerogel->AddMaterial(SiO2, 62.5*perCent);
	//Aerogel->AddMaterial(H2O, 37.4*perCent);
	//Aerogel->AddElement(C, 0.1*perCent);
	
	//Aqui defino las dimensiones de la placa de acero
	G4double plateLength = 1.0 * m;
	G4double plateWidth = 1.0 * m;
   	G4double plateDepth = 10.0 * cm;
	
	G4double energy[2]={1.239841939*eV/0.9, 1.239841939*eV/0.2};
	G4double rindexAerogel[2] = {1.1, 1.1};
	G4double rindexWorld[2]={1.0, 1.0};
	
	
	
	G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");
	
	G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
	mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
	
	worldMat->SetMaterialPropertiesTable(mptWorld);
	
	
////////////////////////////////////////////////////////////////////


//Contruccion del volumen mundo 


	G4Box *solidWorld = new G4Box("solidWorld", 20.0*m, 20.0*m, 20.0*m);
	
	G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	
	G4VisAttributes* invis = new G4VisAttributes();
	 invis->SetVisibility(false);
	logicWorld->SetVisAttributes(invis);
	
	G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);
	
/////////////////////////////////////////////////////////////////

/////////////////Muon Absorber///////////////////



G4Tubs* cil = new G4Tubs("cil", 2.2*m, 2.9*m, 6.6*m, 0.0, 360.0*degree);

G4LogicalVolume *logiccil = new G4LogicalVolume(cil, worldMat, "logiccil");
G4RotationMatrix* rotationc = new G4RotationMatrix();
    rotationc->rotateY(90 * degree);
	
 G4VisAttributes *vis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.5));
	vis->SetVisibility(true);
	logiccil->SetVisAttributes(vis);
	
	G4VPhysicalVolume *physcil = new G4PVPlacement(rotationc, G4ThreeVector(6.0025*m, 0., 0.), logiccil, "physcil", logicWorld, false, 0, true);


/////////////////////////////////////////////////////////


///////////////RPC anillo en ciclo
//Parametros para ciclos///////////////////
//Apotema
///////////////////////////////////////////





//////////////////////////Primer RPC///////////////////////////////




/////////////////Ciclo para Strips///////////////
//Estos strips tienen un ancho de 5 cm y una separacion de 5mm
//separacion entre los centros de los strips
 double separacion = 0.055*m;
//distancia de strips
	double astr1 = 2.91*m ;
	G4Box* solidStripR1 = new G4Box("solidStripR1", 0.025*m, 0.575*m, 0.00005*m);

        // Crear un volumen lógico con el sólido y el material
        logicStripR1 = new G4LogicalVolume(solidStripR1, cobre, "logicStripR1");
for(int x1=0; x1 < 11; x1++) {
 // Crear los volúmenes en un ciclo for
    for (int k = 0; k < 16; k++) {
    double angulogrados = 0+k*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    // l es para las strips que seran 18
    	for (int l = 0; l < 19; l++) {
        // Crear un sólido en forma de caja
         // Calcular las posiciones de los volúmenes
        G4double xPos = x1*1.2*m+0.55*m-l*separacion;  // Separación de 2 veces el tamaño del volumen en el eje x
     
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       
       G4VPhysicalVolume *physStripR1 = new G4PVPlacement(rotation, G4ThreeVector(xPos, astr1*std::sin(anguloradian) , astr1*std::cos(anguloradian)), logicStripR1, "physStripR1", logicWorld, false, 1+l+k*19+304*x1, true);
       	}
    	}
    }


//////////////////Ciclo para primeros acrilicos/////////////
//Distancia del primer acrilico
	double aacrR1 = 2.91055*m ;
	 // Crear un sólido en forma de caja
        G4Box *solidacrilR1 = new G4Box("solidacrilR1", 0.575*m, 0.575*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicacrilR1 = new G4LogicalVolume(solidacrilR1, acrylic, "logicalacrilR1");
 // Crear los volúmenes en un ciclo for
    for (int p2 = 0; p2 < 16; p2++) {
    double angulogrados = 0+p2*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
////ciclo para que se desplace en x 
	for(int x1=0; x1 < 11; x1++) { 
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physacrilR1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, aacrR1*std::sin(anguloradian) , aacrR1*std::cos(anguloradian)), logicacrilR1, "physacrilR1", logicWorld, false, 0, true);
    	   }
       }
    




///////////////////////////////////////////////

///Ciclo para primer pintura//////////////////

//Distancia de la primer pintura
	double apinR1 = 2.9111*m ;
// Crear un sólido en forma de caja
        G4Box *solidpinturaR1 = new G4Box("solidpinturaR1", 0.575*m, 0.575*m, 0.00005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicpinturaR1 = new G4LogicalVolume(solidpinturaR1, grafito, "logicalpinturaR1");
 // Crear los volúmenes en un ciclo for
    for (int p3 = 0; p3 < 16; p3++) {
    double angulogrados = 0+p3*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
    	for(int x1=0; x1 < 11; x1++) {
        // Calcular las posiciones de los volúmenes
        //G4double yPos = 0.415*m-j*separaciony;  // Separación de 2 veces el tamaño del volumen en el eje x
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physpinturaR1 = new G4PVPlacement(rotation, G4ThreeVector(1.2*x1*m, apinR1*std::sin(anguloradian) , apinR1*std::cos(anguloradian)), logicpinturaR1, "physpinturaR1", logicWorld, false, 0, true);
     }
    }


///////////////////////////////////////////


/////Ciclo para primeros vidrios/////////////

//Distancia del primer vidrio
	double avid1R1 = 2.91215*m ;
  // Crear un sólido en forma de caja
        G4Box *solidvid1R1 = new G4Box("solidvid1R1", 0.575*m, 0.575*m, 0.001*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicvid1R1 = new G4LogicalVolume(solidvid1R1, glass, "logicalvid1R1");
 // Crear los volúmenes en un ciclo for
    for (int p4 = 0; p4 < 16; p4++) {
    double angulogrados = 0+p4*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
	for(int x1=0; x1 < 11; x1++) {
        
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physvid1R1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, avid1R1*std::sin(anguloradian) , avid1R1*std::cos(anguloradian)), logicvid1R1, "physvid1R1", logicWorld, false, 0, true);
    }
}


///////////////////////////////////////////



//////////Ciclo para primeros gaps///////////////////


//Distancia del primer gap
	double agap1R1 = 2.91365*m ;
 // Crear un sólido en forma de caja
        G4Box *solidgap1R1 = new G4Box("solidgap1R1", 0.575*m, 0.575*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicgap1R1 = new G4LogicalVolume(solidgap1R1, freonMixture, "logicalgap1R1");
        
 // Crear los volúmenes en un ciclo for
    for (int p5 = 0; p5 < 16; p5++) {
    double angulogrados = 0+p5*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
     for(int x1=0; x1 < 11; x1++) {
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physgap1R1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, agap1R1*std::sin(anguloradian) , agap1R1*std::cos(anguloradian)), logicgap1R1, "physgap1R1", logicWorld, false, 0, true);
    }
}




/////////////////////////////////////////////////////


////////Ciclo para segundo vidrio///////////////

//Distancia del segundo vidrio
	double avid2R1 = 2.91515*m ;
// Crear un sólido en forma de caja
        G4Box *solidvid2R1 = new G4Box("solidvid2R1", 0.575*m, 0.575*m, 0.001*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicvid2R1 = new G4LogicalVolume(solidvid2R1, glass, "logicalvid2R1");
        
 // Crear los volúmenes en un ciclo for
    for (int p6 = 0; p6 < 16; p6++) {
    double angulogrados = 0+p6*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
	for(int x1=0; x1 < 11; x1++) {
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physvid2R1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, avid2R1*std::sin(anguloradian) , avid2R1*std::cos(anguloradian)), logicvid2R1, "physvid2R1", logicWorld, false, 0, true);
    }
}



////////////////////////////////////////////////


/////////////Ciclo para  segundo Gap///////////////

//Distancia del segundo gap
	double agap2R1 = 2.91665*m ;
// Crear un sólido en forma de caja
        G4Box *solidgap2R1 = new G4Box("solidgap2R1", 0.575*m, 0.575*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicgap2R1 = new G4LogicalVolume(solidgap2R1, freonMixture, "logicalgap2R1");
 // Crear los volúmenes en un ciclo for
    for (int p7 = 0; p7 < 16; p7++) {
    double angulogrados = 0+p7*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
     for(int x1=0; x1 < 11; x1++) {
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physgap2R1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, agap2R1*std::sin(anguloradian) , agap2R1*std::cos(anguloradian)), logicgap2R1, "physgap2R1", logicWorld, false, 0, true);
    }
}


////////////////////////////////////////////////////



/////Ciclo para tercer vidrio/////////////////////

//Distancia del tercer vidrio
	double avid3R1 = 2.91815*m ;
// Crear un sólido en forma de caja
        G4Box *solidvid3R1 = new G4Box("solidvid3R1", 0.575*m, 0.575*m, 0.001*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicvid3R1 = new G4LogicalVolume(solidvid3R1, glass, "logicalvid3R1");
        
 // Crear los volúmenes en un ciclo for
    for (int p8 = 0; p8 < 16; p8++) {
    double angulogrados = 0+p8*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
	for(int x1=0; x1 < 11; x1++) {
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physvid3R1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, avid3R1*std::sin(anguloradian) , avid3R1*std::cos(anguloradian)), logicvid3R1, "physvid3R1", logicWorld, false, 0, true);
    }
}


///////////////////////////////////////////////////////////



///////Ciclo para segunda pintura////////////////////////////////


//Distancia de la primer pintura
	double apin2R1 = 2.9192*m ;
// Crear un sólido en forma de caja
        G4Box *solidpintura2R1 = new G4Box("solidpintura2R1", 0.575*m, 0.575*m, 0.00005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicpintura2R1 = new G4LogicalVolume(solidpintura2R1, grafito, "logicalpintura2R1");
 // Crear los volúmenes en un ciclo for
    for (int p9 = 0; p9 < 16; p9++) {
    double angulogrados = 0+p9*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
	for(int x1=0; x1 < 11; x1++) {           
        // Calcular las posiciones de los volúmenes
        //G4double yPos = 0.415*m-j*separaciony;  // Separación de 2 veces el tamaño del volumen en el eje x
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physpintura2R1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, apin2R1*std::sin(anguloradian) , apin2R1*std::cos(anguloradian)), logicpintura2R1, "physpintura2R1", logicWorld, false, 0, true);
    }
}



////////////////////////////////////////////////////////


////////Ciclo para segunda acrilico//////////////////////////////////

//Distancia del primer acrilico
	double aacr2R1 = 2.91975*m ;
// Crear un sólido en forma de caja
        G4Box *solidacril2R1 = new G4Box("solidacril2R1", 0.575*m, 0.575*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicacril2R1 = new G4LogicalVolume(solidacril2R1, acrylic, "logicalacril2R1");
 // Crear los volúmenes en un ciclo for
    for (int p10 = 0; p10 < 16; p10++) {
    double angulogrados = 0+p10*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
	for(int x1=0; x1 < 11; x1++) {  
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physacril2R1 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, aacr2R1*std::sin(anguloradian) , aacr2R1*std::cos(anguloradian)), logicacril2R1, "physacril2R1", logicWorld, false, 0, true);
    }
 }




///////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////
/////////////////////////////////////////////////




//////////////////////////Segundo RPC///////////////////////////////




/////////////////Ciclo para Strips///////////////
// Volumenes en ciclo

/////////////////Ciclo para Strips///////////////
//distancia de strips
	double separaciony = 0.055*m;
	double astr2 = 3.01*m ;
	G4Box* solidStripR2 = new G4Box("solidStripR2", 0.5975*m, 0.025*m, 0.00005*m);

        // Crear un volumen lógico con el sólido y el material
         G4LogicalVolume* logicStripR2 = new G4LogicalVolume(solidStripR2, cobre, "logicStripR2");
//	double l = 2.8038*m;
// Crear los volúmenes en un ciclo for
for(int x1=0; x1 < 11; x1++){//x1 es para ciclar los strips a lo largo del MUONID

    for (int s2 = 0; s2 < 16; s2++) {
    double angulogrados = s2*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    double angulogamarad = 78.75 * M_PI / 180.0;
    	for (int n = 0; n < 21; n++) {
    	//Vamos a usar ley de cosenos
    	 G4double a = 3.0689*m;//no es el apotema, es el lado del triangulo que va del centro al centro de la RPC
    	 G4double b = 0.0262*m+separaciony*n;//esta es la distancia del primer strip al lado que va al centro del rpc 
    	 G4double cosg = std::cos(angulogamarad);//es el coseno del angulo interno del triangulo 
    	 G4double l = sqrt(a*a+b*b-2*a*b*cosg);//ley de cosenos
    	 G4double beta = acos((a*a+l*l-b*b)/(2*a*l));
    	 G4double alpha = (11.25+s2*22.5)* M_PI / 180.0-beta;
        // Crear un sólido en forma de caja
         // Calcular las posiciones de los volúmenes
        G4double yPos = l*std::sin(alpha); 
        G4double zPos = l*std::cos(alpha); 
    G4RotationMatrix* rotationy = new G4RotationMatrix();
    rotationy->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       
       G4VPhysicalVolume *physStripR2 = new G4PVPlacement(rotationy, G4ThreeVector(x1*1.2*m, yPos , zPos), logicStripR2, "physStripR2", logicWorld, false, 1+n+s2*21+x1*336, true);
       	}
    	}
    }


//////////////////Ciclo para primeros acrilicos/////////////
//Distancia del primer acrilico
 // Crear un sólido en forma de caja
        G4Box *solidacrilR2 = new G4Box("solidacrilR2", 0.5975*m, 0.5975*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicacrilR2 = new G4LogicalVolume(solidacrilR2, acrylic, "logicalacrilR2");
	double aacrR2 = 3.01055*m ;
 // Crear los volúmenes en un ciclo for
    for (int a2 = 0; a2 < 16; a2++) {
    double angulogrados = 0+a2*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    for(int x1=0; x1 < 11; x1++) {
        
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physacrilR2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, aacrR2*std::sin(anguloradian) , aacrR2*std::cos(anguloradian)), logicacrilR2, "physacrilR2", logicWorld, false, 0, true);
    	}
    }




///////////////////////////////////////////////

///Ciclo para primer pintura//////////////////

//Distancia de la primer pintura
	double apinR2 = 3.0111*m ;
	// Crear un sólido en forma de caja
        G4Box *solidpinturaR2 = new G4Box("solidpinturaR2", 0.575*m, 0.575*m, 0.00005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicpinturaR2 = new G4LogicalVolume(solidpinturaR2, grafito, "logicalpinturaR2");
 // Crear los volúmenes en un ciclo for
    for (int a3 = 0; a3 < 16; a3++) {
    double angulogrados = 0+a3*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    for(int x1=0; x1 < 11; x1++) {
    
        // Crear un sólido en forma de caja
        G4Box *solidpinturaR2 = new G4Box("solidpinturaR2", 0.575*m, 0.575*m, 0.00005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicpinturaR2 = new G4LogicalVolume(solidpinturaR2, grafito, "logicalpinturaR2");
       
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physpinturaR2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, apinR2*std::sin(anguloradian) , apinR2*std::cos(anguloradian)), logicpinturaR2, "physpinturaR2", logicWorld, false, 0, true);
     }
    }


///////////////////////////////////////////


/////Ciclo para primeros vidrios/////////////

//Distancia del primer vidrio
	double avid1R2 = 3.01215*m ;
	// Crear un sólido en forma de caja
        G4Box *solidvid1R2 = new G4Box("solidvid1R2", 0.5975*m, 0.5975*m, 0.001*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicvid1R2 = new G4LogicalVolume(solidvid1R2, glass, "logicalvid1R2");
        
 // Crear los volúmenes en un ciclo for
    for (int a4 = 0; a4 < 16; a4++) {
    double angulogrados = 0+a4*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
    for(int x1=0; x1 < 11; x1++) {
      
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physvid1R2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, avid1R2*std::sin(anguloradian) , avid1R2*std::cos(anguloradian)), logicvid1R2, "physvid1R2", logicWorld, false, 0, true);
    }
}


///////////////////////////////////////////



//////////Ciclo para primeros gaps///////////////////


//Distancia del primer gap
	double agap1R2 = 3.01365*m ;
	 // Crear un sólido en forma de caja
        G4Box *solidgap1R2 = new G4Box("solidgap1R2", 0.5975*m, 0.5975*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicgap1R2 = new G4LogicalVolume(solidgap1R2, freonMixture, "logicalgap1R2");
 // Crear los volúmenes en un ciclo for
    for (int a5 = 0; a5 < 16; a5++) {
    double angulogrados = 0+a5*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
    for(int x1=0; x1 < 11; x1++) {
        
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physgap1R2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, agap1R2*std::sin(anguloradian) , agap1R2*std::cos(anguloradian)), logicgap1R2, "physgap1R2", logicWorld, false, 0, true);
    }
   }





/////////////////////////////////////////////////////


////////Ciclo para segundo vidrio///////////////

//Distancia del segundo vidrio
	double avid2R2 = 3.01515*m ;
	// Crear un sólido en forma de caja
        G4Box *solidvid2R2 = new G4Box("solidvid2R2", 0.5975*m, 0.5975*m, 0.001*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicvid2R2 = new G4LogicalVolume(solidvid2R2, glass, "logicalvid2R2");
 // Crear los volúmenes en un ciclo for
    for (int a6 = 0; a6 < 16; a6++) {
    double angulogrados = 0+a6*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
    for(int x1=0; x1 < 11; x1++) {
        
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physvid2R2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, avid2R2*std::sin(anguloradian) , avid2R2*std::cos(anguloradian)), logicvid2R2, "physvid2R2", logicWorld, false, 0, true);
    }
}



////////////////////////////////////////////////


/////////////Ciclo para  segundo Gap///////////////

//Distancia del segundo gap
	double agap2R2 = 3.01665*m ;
// Crear un sólido en forma de caja
        G4Box *solidgap2R2 = new G4Box("solidgap2R2", 0.5975*m, 0.5975*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicgap2R2 = new G4LogicalVolume(solidgap2R2, freonMixture, "logicalgap2R2");
 // Crear los volúmenes en un ciclo for
    for (int a7 = 0; a7 < 16; a7++) {
    double angulogrados = 0+a7*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
    for(int x1=0; x1 < 11; x1++) {
        
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physgap2R2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, agap2R2*std::sin(anguloradian) , agap2R2*std::cos(anguloradian)), logicgap2R2, "physgap2R2", logicWorld, false, 0, true);
    }
 }


////////////////////////////////////////////////////



/////Ciclo para tercer vidrio/////////////////////

//Distancia del tercer vidrio
	double avid3R2 = 3.01815*m ;
	// Crear un sólido en forma de caja
        G4Box *solidvid3R2 = new G4Box("solidvid3R2", 0.5975*m, 0.5975*m, 0.001*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicvid3R2 = new G4LogicalVolume(solidvid3R2, glass, "logicalvid3R2");
 // Crear los volúmenes en un ciclo for
    for (int a8 = 0; a8 < 16; a8++) {
    double angulogrados = 0+a8*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    
    for(int x1=0; x1 < 11; x1++) {
        
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physvid3R2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, avid3R2*std::sin(anguloradian) , avid3R2*std::cos(anguloradian)), logicvid3R2, "physvid3R2", logicWorld, false, 0, true);
    }
}




///////////////////////////////////////////////////////////



///////Ciclo para segunda pintura////////////////////////////////


//Distancia de la primer pintura
	double apin2R2 = 3.0192*m ;
	 // Crear un sólido en forma de caja
        G4Box *solidpintura2R2 = new G4Box("solidpintura2R2", 0.5975*m, 0.5975*m, 0.00005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicpintura2R2 = new G4LogicalVolume(solidpintura2R2, grafito, "logicalpintura2R2");
 // Crear los volúmenes en un ciclo for
    for (int a9 = 0; a9 < 16; a9++) {
    double angulogrados = 0+a9*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    for(int x1=0; x1 < 11; x1++) {
    
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physpintura2R2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, apin2R2*std::sin(anguloradian) , apin2R2*std::cos(anguloradian)), logicpintura2R2, "physpintura2R2", logicWorld, false, 0, true);
    }
}



////////////////////////////////////////////////////////


////////Ciclo para segunda acrilico//////////////////////////////////

//Distancia del primer acrilico
	double aacr2R2 = 3.01975*m ;
        // Crear un sólido en forma de caja
        G4Box *solidacril2R2 = new G4Box("solidacril2R2", 0.5975*m, 0.5975*m, 0.0005*m);

        // Crear un volumen lógico con el sólido y el material
        G4LogicalVolume *logicacril2R2 = new G4LogicalVolume(solidacril2R2, acrylic, "logicalacril2R2");
 // Crear los volúmenes en un ciclo for
    for (int a10 = 0; a10 < 16; a10++) {
    double angulogrados = 0+a10*22.5;
    // Convertir el ángulo de grados a radianes
    double anguloradian = angulogrados * M_PI / 180.0;
    for(int x1=0; x1 < 11; x1++) {
    
        
// Crear una matriz de rotación para la rotación deseada
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(angulogrados * degree);
        // Colocar el volumen en la posición calculada
       G4VPhysicalVolume *physacril2R2 = new G4PVPlacement(rotation, G4ThreeVector(x1*1.2*m, aacr2R2*std::sin(anguloradian) , aacr2R2*std::cos(anguloradian)), logicacril2R2, "physacril2R2", logicWorld, false, 0, true);
    }
}




///////////////////////////////////////////////////////////////

 
	return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
	SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");
	
	logicStripR1->SetSensitiveDetector(sensDet);
}
