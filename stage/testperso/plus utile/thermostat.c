//Declaration des etats et initialisation
unsigned int ES,EP;	// Variables Etat suivant et Etat present
ES=EP=0;

/* Correspondance entre les nom d'etats et le vecteur entier : 
*  E0 : 0
*  E1 : 1
*/
// Debut du bloc F de la realisation


switch(EP) {
	// code pour Ep = E0
	case 0 : 
		if(cdt1)
		{
			ES = 1;
		}
		break;
	// code pour Ep = E1
	case 1 : 
		if(cdt2)
		{
			ES = 0;
		}
		break;
}
EP=ES;
//Fin du bloc F de la realisation

//Debut de realisation du bloc G : 

if( EP=0 )
	marche = 1;
else
	marche = 0;

if( EP=1 )
	arret = 1;
else
	arret = 0;


//Fin de realisation du bloc G  
