// =======================================================================
// Purpose:     header for Covid model population data with new kernels of ParMooN
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 04.04.2020
// =======================================================================

// enum STATE {AN, AP, AR, AS, BR, CH,CT,DL,GA,GJ,HR,HP,JK,JH,KA,KL,LA,LD,MP,MH,MN,OR,PY,PB,RJ,SK,TN,TG,TR,UP,UT,WB};

// string enum_to_string(int type) {
//    switch(type) {
//       case 0:
//          return "AN";
//       case 1:
//          return "AP";
//       case 2:
//          return "AR";
//       case 3:
//          return "AS";
//       case 4:
//          return "BR";
//       case 5:
//          return "CH";
//       case 6:
//          return "CT";
//       case 7:
//          return "DL";
//       case 8:
//          return "GA";
//       case 9:
//          return "GJ";
//       case 10:
//          return "HR";
//       case 11:
//          return "HP";
//       case 12:
//          return "JK";
//       case 13:
//          return "JH";
//       case 14:
//          return "KA";
//       case 15:
//          return "KL";
//       case 16:
//          return "LA";         
//       case 17:
//          return "LD";
//       case 18:
//          return "MP";
//       case 19:
//          return "MH";   
//       case 20:
//          return "MN"; 
//       case 21:
//          return "OR";
//       case 22:
//          return "PY";
//       case 23:
//          return "PB";
//       case 24:
//          return "RJ";
//       case 25:
//          return "SK";
//       case 26:
//          return "TN";
//       case 27:
//          return "TG";
//       case 28:
//          return "TR";
//       case 29:
//          return "UP";
//       case 30:
//          return "UT";
//       case 31:
//          return "WB";
//       default:
//          cout <<"Invalid District"<<endl;
//          exit(0);
//    }
// }

// Karnataka
enum STATE {Bagalkote, Ballari, Belagavi, BengaluruRural, BengaluruUrban, Bidar, Chamarajanagara, 
            Chikkaballapura, Chikkamagaluru, Chitradurga, DakshinaKannada, Davanagere, Dharwad, Gadag, 
            Hassan, Haveri, Kalaburagi, Kodagu, Kolar, Koppal, Mandya, Mysuru, Raichur, Ramanagara, 
            Shivamogga, Tumakuru, Udupi, UttaraKannada, Vijayapura, Yadgir, unknown, unknown1};

//   //CIR Wave-I
//   double State_CIR2[32] = {, 1., 1.};
  
string enum_to_string(int type) {
   switch(type) {
case 0:
    return "Bagalkote";
case 1:
    return "Ballari";
case 2:
    return "Belagavi";
case 3:
    return "BengaluruRural";
case 4:
    return "BengaluruUrban";
case 5:
    return "Bidar";
case 6:
    return "Chamarajanagara";
case 7:
    return "Chikkaballapura";
case 8:
    return "Chikkamagaluru";
case 9:
    return "Chitradurga";
case 10:
    return "DakshinaKannada";
case 11:
    return "Davanagere";
case 12:
    return "Dharwad";
case 13:
    return "Gadag";
case 14:
    return "Hassan";
case 15:
    return "Haveri";
case 16:
    return "Kalaburagi";
case 17:
    return "Kodagu";
case 18:
    return "Kolar";
case 19:
    return "Koppal";
case 20:
    return "Mandya";
case 21:
    return "Mysuru";
case 22:
    return "Raichur";
case 23:
    return "Ramanagara";
case 24:
    return "Shivamogga";
case 25:
    return "Tumakuru";
case 26:
    return "Udupi";
case 27:
    return "UttaraKannada";
case 28:
    return "Vijayapura";
case 29:
    return "Yadgir";
case 30:
    return "unknown";
case 31:
    return "unknown1";
      default:
         cout <<"Invalid District"<<endl;
         exit(0);
   }
}

// // Jarkhand Distric data
// enum STATE {Bokaro, Chatra, Deoghar, Dhanbad, Dumka, EastSinghbhum, Garhwa, Giridih, Godda,
//             Gumla, Hazaribagh, Jamtara, Khunti, Koderma, Latehar, Lohardaga, Pakur, Palamu,
//             Ramgarh, Ranchi, Sahibganj, SaraikelaKharsawan, Simdega, WestSinghbhum, unknown1,
//             unknown2, unknown3, unknown4, unknown5, unknown6, unknown7, unknown8};

        
// string enum_to_string(int type) {
//    switch(type) {
// case 0:
//     return "Bokaro";
// case 1:
//     return "Chatra";
// case 2:
//     return "Deoghar";
// case 3:
//     return "Dhanbad";
// case 4:
//     return "Dumka";
// case 5:
//     return "EastSinghbhum";
// case 6:
//     return "Garhwa";
// case 7:
//     return "Giridih";
// case 8:
//     return "Godda";
// case 9:
//     return "Gumla";
// case 10:
//     return "Hazaribagh";
// case 11:
//     return "Jamtara";
// case 12:
//     return "Khunti";
// case 13:
//     return "Koderma";
// case 14:
//     return "Latehar";
// case 15:
//     return "Lohardaga";
// case 16:
//     return "Pakur";
// case 17:
//     return "Palamu";
// case 18:
//     return "Ramgarh";
// case 19:
//     return "Ranchi";
// case 20:
//     return "Sahibganj";
// case 21:
//     return "SaraikelaKharsawan";
// case 22:
//     return "Simdega";
// case 23:
//     return "WestSinghbhum";
// case 24:
//     return "unknown1";
// case 25:
//     return "unknown2";
// case 26:
//     return "unknown3";
// case 27:
//     return "unknown4";
// case 28:
//     return "unknown5";
// case 29:
//     return "unknown6";
// case 30:
//     return "unknown7";
// case 31:
//     return "unknown8";
//       default:
//          cout <<"Invalid District"<<endl;
//          exit(0);
//    }
// }    

// Tamilnadu
// enum STATE {Ariyalur, Chengalpattu, Chennai, Coimbatore, Cuddalore, Dharmapuri, Dindigul, Erode, Kallakurichi,
//             Kancheepuram, Kanyakumari, Karur, Krishnagiri, Madurai, Nagapattinam, Namakkal, Nilgiris, Perambalur, Pudukkottai,
//             Ramanathapuram, Ranipet, Salem, Sivaganga, Tenkasi, Thanjavur, Theni, Thiruvallur, Thiruvarur, Thoothukkudi, 
//             Tiruchirappalli, Tirunelveli, Tirupathur, Tiruppur, Tiruvannamalai, Vellore, Viluppuram, Virudhunagar};

// string enum_to_string(int type) {
//    switch(type) {
// case 0:
//     return "Ariyalur";
// case 1:
//     return "Chengalpattu";
// case 2:
//     return "Chennai";
// case 3:
//     return "Coimbatore";
// case 4:
//     return "Cuddalore";
// case 5:
//     return "Dharmapuri";
// case 6:
//     return "Dindigul";
// case 7:
//     return "Erode";
// case 8:
//     return "Kallakurichi";
// case 9:
//     return "Kancheepuram";
// case 10:
//     return "Kanyakumari";
// case 11:
//     return "Karur";
// case 12:
//     return "Krishnagiri";
// case 13:
//     return "Madurai";
// case 14:
//     return "Nagapattinam";
// case 15:
//     return "Namakkal";
// case 16:
//     return "Nilgiris";
// case 17:
//     return "Perambalur";
// case 18:
//     return "Pudukkottai";
// case 19:
//     return "Ramanathapuram";
// case 20:
//     return "Ranipet";
// case 21:
//     return "Salem";
// case 22:
//     return "Sivaganga";
// case 23:
//     return "Tenkasi";
// case 24:
//     return "Thanjavur";
// case 25:
//     return "Theni";
// case 26:
//     return "Thiruvallur";
// case 27:
//     return "Thiruvarur";
// case 28:
//     return "Thoothukkudi";
// case 29:
//     return "Tiruchirappalli";
// case 30:
//     return "Tirunelveli";
// case 31:
//     return "Tirupathur";
// case 32:
//     return "Tiruppur";
// case 33:
//     return "Tiruvannamalai";
// case 34:
//     return "Vellore";
// case 35:
//     return "Viluppuram";
// case 36:
//     return "Virudhunagar";
// case 37:
//     return "unknown1"; 
// case 38:
//     return "unknown2";
// case 39:
//     return "unknown3";                         
//       default:
//          cout << type <<" : Invalid District"<<endl;
//          exit(0);
//    }
// }
 
// // Pondichery
// enum STATE {Karaikal, Mahe, Puducherry, Yanam};

// string enum_to_string(int type) {
//    switch(type) {
// case 0:
//     return "Karaikal";
// case 1:
//     return "Mahe";
// case 2:
//     return "Puducherry";
// case 3:
//     return "Yanam";
// default:
//          cout << type <<" : Invalid District"<<endl;
//          exit(0);
//    }
// }
